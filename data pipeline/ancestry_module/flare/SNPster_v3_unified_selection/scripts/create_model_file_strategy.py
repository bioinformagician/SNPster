#!/usr/bin/env python3
"""Create an initial FLARE2 model from copying-profile panel probabilities.

This is a strategy-search extension of the Browning Lab FLARE create_model_file
idea.  It deliberately keeps the FLARE model-file output common while allowing
multiple clustering strategies to compete:

  * gaussian_mixture
  * bayesian_gaussian_mixture
  * kmeans
  * minibatch_kmeans

Any strategy may optionally cluster in PCA space.  Cluster centres are mapped
back to the original copying-probability simplex before the FLARE model is
written.  FLARE then performs the same downstream model refinement for every
strategy, so the benchmark comparison remains apples-to-apples.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.mixture import BayesianGaussianMixture, GaussianMixture


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Expanded model JSON produced by generate_model_configs.py")
    ap.add_argument("--panels", required=True)
    ap.add_argument("--output-prefix", required=True)
    return ap.parse_args()


def read_panels(path: str):
    with open(path) as panelfile:
        nrefhaps = int(panelfile.readline().split("=")[1])
        refpops = panelfile.readline().split()[4:]
        npanel = len(refpops)

        panelprobs = []
        rhovals = []
        panelprobschrom = []
        rhovalschrom = []
        oldlocindex = -1
        hapindices = []
        oldchrom = None
        oldhapindices = []

        for line in panelfile:
            if not line or line[0] == "#":
                continue
            bits = line.split()
            locindex = int(bits[0])
            chrom = bits[1].split(":")[0]
            hapindex = int(bits[2])
            thisrho = float(bits[3])

            if chrom != oldchrom:
                oldchrom = chrom
                if panelprobschrom:
                    panelprobs.extend(panelprobschrom)
                    rhovals.extend(rhovalschrom)
                panelprobschrom = []
                rhovalschrom = []
                oldlocindex = -1

            if locindex != oldlocindex:
                oldlocindex = locindex
                panelprobschrom.append([])
                rhovalschrom.append([])
                if hapindices and oldhapindices and any(x != y for x, y in zip(hapindices, oldhapindices)):
                    raise SystemExit("haplotype indices inconsistent across positions")
                oldhapindices = list(hapindices)
                hapindices = []

            hapindices.append(hapindex)
            probs = [float(x) for x in bits[4:]]
            if len(probs) != npanel:
                raise SystemExit("number of reference panels varies in .panels file")
            panelprobschrom[-1].append(probs)
            rhovalschrom[-1].append(thisrho)

        panelprobs.extend(panelprobschrom)
        rhovals.extend(rhovalschrom)

    if not panelprobs:
        raise SystemExit("No copying profiles found in panels file")

    nhaps = len(panelprobs[0])
    if any(len(loc) != nhaps for loc in panelprobs):
        raise SystemExit("number of target haplotypes varies across locations")

    X = np.asarray([hap for loc in panelprobs for hap in loc], dtype=np.float64)
    rho_vec = np.asarray([hap for loc in rhovals for hap in loc], dtype=np.float64)
    return nrefhaps, refpops, len(panelprobs), nhaps, X, rho_vec


def optional(cfg: dict[str, Any], key: str, default=None):
    value = cfg.get(key, default)
    return default if value is None else value


def make_clusterer(cfg: dict[str, Any], n_components: int):
    family = cfg["family"]
    seed = int(cfg.get("model_seed", 0))
    n_init = int(cfg.get("n_init", 1))
    max_iter = int(cfg.get("max_iter", 300))
    tol = float(cfg.get("tol", 1e-3))

    if family == "gaussian_mixture":
        return GaussianMixture(
            n_components=n_components,
            covariance_type=cfg.get("covariance_type", "diag"),
            init_params=cfg.get("init_params", "k-means++"),
            n_init=n_init,
            max_iter=max_iter,
            tol=tol,
            reg_covar=float(cfg.get("reg_covar", 1e-6)),
            random_state=seed,
        )

    if family == "bayesian_gaussian_mixture":
        return BayesianGaussianMixture(
            n_components=n_components,
            covariance_type=cfg.get("covariance_type", "diag"),
            init_params=cfg.get("init_params", "k-means++"),
            n_init=n_init,
            max_iter=max_iter,
            tol=tol,
            reg_covar=float(cfg.get("reg_covar", 1e-6)),
            weight_concentration_prior_type=cfg.get(
                "weight_concentration_prior_type", "dirichlet_process"
            ),
            weight_concentration_prior=cfg.get("weight_concentration_prior", None),
            mean_precision_prior=cfg.get("mean_precision_prior", None),
            random_state=seed,
        )

    if family == "kmeans":
        return KMeans(
            n_clusters=n_components,
            init=cfg.get("kmeans_init", "k-means++"),
            n_init=n_init,
            max_iter=max_iter,
            tol=tol,
            random_state=seed,
            algorithm=cfg.get("kmeans_algorithm", "lloyd"),
        )

    if family == "minibatch_kmeans":
        return MiniBatchKMeans(
            n_clusters=n_components,
            init=cfg.get("kmeans_init", "k-means++"),
            n_init=n_init,
            max_iter=max_iter,
            tol=tol,
            batch_size=int(cfg.get("batch_size", 4096)),
            reassignment_ratio=float(cfg.get("reassignment_ratio", 0.01)),
            random_state=seed,
        )

    raise SystemExit(f"Unsupported model family: {family}")


def cluster_centers(model):
    if hasattr(model, "means_"):
        return np.asarray(model.means_, dtype=np.float64)
    if hasattr(model, "cluster_centers_"):
        return np.asarray(model.cluster_centers_, dtype=np.float64)
    raise SystemExit(f"Cannot obtain centres from {type(model).__name__}")


def model_weights(model, k: int):
    if hasattr(model, "weights_"):
        return np.asarray(model.weights_, dtype=np.float64)
    return np.full(k, 1.0 / k, dtype=np.float64)


def compact_components(cfg, labels, centers, weights, cluster_X):
    requested_k = len(centers)
    counts = np.bincount(labels, minlength=requested_k)
    compact = bool(cfg.get("compact_inactive", False))
    min_weight = float(cfg.get("min_component_weight", 0.0))

    if compact:
        active = (counts > 0) & (weights >= min_weight)
        # Never discard a component that has meaningful assigned mass simply
        # because the variational weight is microscopically below threshold.
        frac = counts / max(1, counts.sum())
        active |= frac >= max(min_weight, 1e-8)
    else:
        active = counts > 0

    active_idx = np.where(active)[0]
    if len(active_idx) < 2:
        raise SystemExit(
            f"Only {len(active_idx)} active components remain from requested K={requested_k}"
        )

    if len(active_idx) != requested_k:
        # Reassign observations belonging to discarded tiny/empty Bayesian
        # components to the nearest retained centre in clustering space.
        inactive_mask = ~np.isin(labels, active_idx)
        if inactive_mask.any():
            pts = cluster_X[inactive_mask]
            keep_centers = centers[active_idx]
            # Small K, so the explicit distance matrix is acceptable here.
            d2 = ((pts[:, None, :] - keep_centers[None, :, :]) ** 2).sum(axis=2)
            labels = labels.copy()
            labels[inactive_mask] = active_idx[np.argmin(d2, axis=1)]

        remap = {old: new for new, old in enumerate(active_idx)}
        labels = np.asarray([remap[int(x)] for x in labels], dtype=np.int32)
        centers = centers[active_idx]
        weights = weights[active_idx]
        weights = weights / weights.sum()
        counts = np.bincount(labels, minlength=len(active_idx))
    else:
        labels = labels.astype(np.int32, copy=False)

    return labels, centers, weights, counts, active_idx


def main():
    args = parse_args()
    cfg = json.loads(Path(args.config).read_text())
    prefix = Path(args.output_prefix)

    nrefhaps, refpops, nloc, nhaps, X, rho_vec = read_panels(args.panels)
    npanel = len(refpops)
    requested_k = int(cfg["k"])

    # FLARE's original clustering construction drops one simplex coordinate to
    # remove the exact sum-to-one dependence before clustering.
    whichdrop = int(np.argmax(X.sum(axis=0)))
    transformed = np.delete(X, whichdrop, axis=1)

    pca_components = int(cfg.get("pca_components", 0) or 0)
    pca = None
    cluster_X = transformed
    if pca_components > 0:
        max_pc = min(transformed.shape[1], transformed.shape[0] - 1)
        if pca_components > max_pc:
            raise SystemExit(
                f"pca_components={pca_components} exceeds maximum possible {max_pc}"
            )
        pca = PCA(
            n_components=pca_components,
            svd_solver="randomized",
            whiten=bool(cfg.get("pca_whiten", False)),
            random_state=int(cfg.get("model_seed", 0)),
        )
        cluster_X = pca.fit_transform(transformed)

    model = make_clusterer(cfg, requested_k)
    labels = model.fit_predict(cluster_X) if hasattr(model, "fit_predict") else None
    if labels is None:
        model.fit(cluster_X)
        labels = model.predict(cluster_X)

    centers_cluster = cluster_centers(model)
    weights = model_weights(model, requested_k)
    labels, centers_cluster, weights, counts, active_idx = compact_components(
        cfg, np.asarray(labels), centers_cluster, weights, cluster_X
    )
    effective_k = len(centers_cluster)

    # Map centres back out of PCA (if used), then back onto the complete
    # copying-probability simplex by restoring the coordinate FLARE drops.
    centers_reduced = pca.inverse_transform(centers_cluster) if pca is not None else centers_cluster
    extra = 1.0 - centers_reduced.sum(axis=1)
    param = np.insert(centers_reduced, whichdrop, extra, axis=1)
    param = np.maximum(param, 0.0)
    row_sums = param.sum(axis=1, keepdims=True)
    if np.any(row_sums <= 0):
        raise SystemExit("A cluster centre collapsed to a zero vector after simplex projection")
    param = param / row_sums

    rho = []
    for i in range(effective_k):
        vals = rho_vec[labels == i]
        if len(vals) == 0:
            raise SystemExit(f"Effective component {i} has no assigned observations")
        rho.append(float(np.mean(vals)))

    T = 10
    model_path = Path(str(prefix) + ".model")
    with model_path.open("w") as outfile:
        print("# configurable FLARE2 strategy model", file=outfile)
        print("# " + json.dumps(cfg, sort_keys=True), file=outfile)
        print("\t".join(f"anc_{i}" for i in range(effective_k)), file=outfile)
        print(file=outfile)
        print("\t".join(refpops), file=outfile)
        print(file=outfile)
        print(T, file=outfile)
        print(file=outfile)
        print("\t".join([str(1.0 / effective_k)] * effective_k), file=outfile)
        print(file=outfile)
        for row in param:
            print("\t".join(map(str, row)), file=outfile)
        print(file=outfile)
        mylambda = 1 / (np.log(nrefhaps) + 0.5)
        theta = mylambda / (2 * mylambda + 2 * nrefhaps)
        for _ in range(effective_k):
            print("\t".join([str(theta)] * npanel), file=outfile)
        print(file=outfile)
        print("\t".join(map(str, rho)), file=outfile)

    summary = {
        "strategy_name": cfg.get("strategy_name"),
        "family": cfg["family"],
        "requested_k": requested_k,
        "effective_k": effective_k,
        "active_original_components": [int(x) for x in active_idx],
        "component_counts": [int(x) for x in counts],
        "component_weights": [float(x) for x in weights],
        "n_locations": nloc,
        "n_haplotypes": nhaps,
        "n_observations": int(X.shape[0]),
        "n_reference_panels": npanel,
        "n_reference_haplotypes": nrefhaps,
        "which_panel_column_dropped": whichdrop,
        "pca_components": pca_components,
        "pca_explained_variance": (
            float(np.sum(pca.explained_variance_ratio_)) if pca is not None else None
        ),
        "converged": bool(getattr(model, "converged_", True)),
        "n_iter": int(getattr(model, "n_iter_", -1)),
        "lower_bound": (
            float(getattr(model, "lower_bound_"))
            if hasattr(model, "lower_bound_") and math.isfinite(float(getattr(model, "lower_bound_")))
            else None
        ),
        "inertia": (
            float(getattr(model, "inertia_")) if hasattr(model, "inertia_") else None
        ),
        "config": cfg,
    }
    Path(str(prefix) + ".json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    Path(str(prefix) + ".log").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    if not summary["converged"]:
        raise SystemExit("Clustering model did not converge")


if __name__ == "__main__":
    main()

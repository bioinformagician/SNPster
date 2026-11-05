#!/usr/bin/env python3
"""Shared donor-purged splits and accuracy reports for phase-independent sweeps."""

import hashlib

import numpy as np
import pandas as pd

import model_family_gnomix as adapter
from report_phasefree_accuracy import accuracy_summary
from screen_model_families import unified_seven
from screen_v1_decoders import purged_folds, training_weights


SOURCES = list(adapter.SOURCES)
FEATURES = [f"{view}_{source}" for view in ("classifier", "frequency") for source in SOURCES]
TRUTH = [f"truth_{source}" for source in SOURCES]


def stable_features(values, decimals=12):
    values = np.asarray(values, float)
    if (values.ndim != 2 or values.shape[1] != len(FEATURES) or not len(values)
            or not np.isfinite(values).all() or np.any(values < -1e-12)
            or np.any(values > 1 + 1e-12)
            or not np.allclose(values[:, :7].sum(axis=1), 1, atol=1e-10, rtol=0)
            or not np.allclose(values[:, 7:].sum(axis=1), 1, atol=1e-10, rtol=0)):
        raise ValueError("Expected finite paired seven-source classifier/frequency features")
    return np.round(values, decimals=decimals)


def donor_purged_splits(metadata, donors, pure_folds, seed):
    if (metadata.sample_id.duplicated().any() or set(metadata.track) != {"pure", "pedigree"}
            or pure_folds < 2):
        raise ValueError("Invalid training-only split metadata")
    pure = metadata[metadata.track.eq("pure")]
    pedigree = metadata[metadata.track.eq("pedigree")]
    if (set(pure.source_group) != set(SOURCES) or set(pedigree.family_id) != set(donors)
            or not set().union(*donors.values()) <= set(pure.sample_id)):
        raise ValueError("Incomplete training source/founder/family metadata")
    splits = []
    for family, train, test, info in purged_folds(pedigree.reset_index(drop=True), donors):
        held = set(info["test_donors"])
        train_ids = set(pedigree.iloc[train].sample_id) | set(pure.loc[~pure.sample_id.isin(held), "sample_id"])
        test_ids = set(pedigree.iloc[test].sample_id)
        splits.append(dict(fold_id=f"family_{family}", track="pedigree",
                           heldout_donors=sorted(held), train_ids=sorted(train_ids),
                           test_ids=sorted(test_ids)))
    assignment = {}
    for _, group in pure.groupby("source_group"):
        if len(group) < pure_folds:
            raise ValueError("Too few pure training donors for source-stratified folds")
        ordered = sorted(group.sample_id, key=lambda identity: hashlib.sha256(
            f"{seed}:{identity}".encode()).digest())
        assignment.update({identity: index % pure_folds for index, identity in enumerate(ordered)})
    for fold in range(pure_folds):
        held = {identity for identity, value in assignment.items() if value == fold}
        safe_families = {family for family, identities in donors.items() if not identities & held}
        train_ids = set(pure.loc[~pure.sample_id.isin(held), "sample_id"])
        train_ids |= set(pedigree.loc[pedigree.family_id.isin(safe_families), "sample_id"])
        splits.append(dict(fold_id=f"pure_{fold}", track="pure", heldout_donors=sorted(held),
                           train_ids=sorted(train_ids), test_ids=sorted(held)))
    seen = []
    for split in splits:
        train = metadata[metadata.sample_id.isin(split["train_ids"])]
        test = metadata[metadata.sample_id.isin(split["test_ids"])]
        held = set(split["heldout_donors"])
        train_donors = set(train.loc[train.track.eq("pure"), "sample_id"])
        train_donors |= set().union(*(donors[family] for family in
                                    train.loc[train.track.eq("pedigree"), "family_id"].unique()))
        if (not len(train) or not len(test) or train_donors & held
                or set(train.sample_id) & set(test.sample_id)
                or set(train.loc[train.track.eq("pure"), "source_group"]) != set(SOURCES)
                or not train.track.eq("pedigree").any()):
            raise ValueError(f"Invalid donor-purged training fold: {split['fold_id']}")
        split.update(n_train=len(train), n_test=len(test), n_train_pure=int(train.track.eq("pure").sum()),
                     n_train_pedigree=int(train.track.eq("pedigree").sum()), donor_overlap=0)
        seen.extend(split["test_ids"])
    if len(seen) != len(metadata) or set(seen) != set(metadata.sample_id):
        raise ValueError("Every clean training person must have exactly one out-of-fold prediction")
    return splits


def fold_weights(metadata):
    pure = metadata[metadata.track.eq("pure")]
    pedigree = metadata[metadata.track.eq("pedigree")]
    if not len(pure) or not len(pedigree) or set(pure.source_group) != set(SOURCES):
        raise ValueError("Fold weighting requires pure-source coverage and pedigree training")
    weights = pd.Series(training_weights(pure, pedigree),
                        index=pd.concat([pure, pedigree]).sample_id)
    return weights.loc[metadata.sample_id].to_numpy(float)


def summarize(metadata, predictions, targets, gates):
    if len(metadata) != len(predictions) or metadata.sample_id.duplicated().any():
        raise ValueError("Validation metadata/prediction axes differ")
    frame = metadata.copy()
    for index, source in enumerate(SOURCES):
        frame[f"pred_{source}"] = predictions[:, index]
    unified = unified_seven(frame, gates)
    result, breakdown = dict(unified=unified), []
    for track, group in frame.groupby("track"):
        values = group[[f"pred_{source}" for source in SOURCES]].to_numpy(float)
        truth = group[TRUTH].to_numpy(float)
        summary, _ = accuracy_summary(truth, values, targets)
        summary["off_target_mass"] = float((values * (truth < gates["absent_ancestry_max"])).sum(axis=1).mean())
        result[track] = summary
        for column in (("source_group",) if track == "pure" else ("scenario",)):
            for name, subset in group.groupby(column):
                summary, _ = accuracy_summary(
                    subset[TRUTH].to_numpy(float),
                    subset[[f"pred_{source}" for source in SOURCES]].to_numpy(float), targets)
                breakdown.append(dict(track=track, grouping=column, value=name, **summary))
    result["breakdown"] = breakdown
    return result


def ranking_key(candidate, result):
    return (result["unified"]["balanced_mean_tv"], result["pedigree"]["mean_tv"],
            -result["pedigree"]["fraction_all_components_within_5pp"],
            -result["pedigree"]["fraction_all_components_within_10pp"], candidate)

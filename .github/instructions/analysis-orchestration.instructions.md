---
description: "Use when planning, implementing or launching SNPster analyses, ancestry benchmarks, model sweeps or validation experiments. Require containerized Nextflow orchestration, graphical pipeline progress, Docker logs and scoped cache cleanup."
applyTo: "data_analysis/**"
---

# Analysis orchestration

- Use Nextflow DSL2 for orchestration and parallelization of new or expanded multi-step analysis workflows, model sweeps and benchmarking experiments in this project.
- Run the Nextflow orchestrator itself in a named Docker container, not just its individual tasks. Provide the actual container name and `docker logs -f <container_name>` command at launch. Send overall progress, failures, verification and cleanup status to its stdout/stderr. Do not auto-remove the controller before its logs and exit status have been preserved.
- For future runs, enable Nextflow's graphical/live process progress display with `-ansi-log true`, allocating a controller TTY (`docker run -t`) as needed. Verify that the live pipeline overview is visible through the supplied monitoring command. Also save a whole-pipeline DAG (`-with-dag pipeline.html`), HTML report and timeline. Retain the plain Nextflow log and task trace alongside the graphical view for durable debugging and auditing. Do not retrofit these display settings into an already-running frozen experiment.
- Pin the controller image and Nextflow version, verify its task-container engine access, and preserve host/container absolute path consistency for task mounts. Persist the dedicated launch/history, work and published-output directories outside the controller so failed runs remain resumable. Verify a real task-container smoke test before launching a full run.
- Keep Python, R and shell scripts as task implementations. Do not introduce another top-level Python process pool or bespoke shell scheduler in place of Nextflow.
- Make reusable preparation, phasing, model execution, aggregation and verification separate cacheable processes where appropriate. Declare inputs, outputs, pinned containers, seeds and CPU/memory requirements explicitly; bound concurrency to the available resource budget.
- Use Nextflow task caching and `-resume` for reproducible reuse. Preserve scientific input hashes, donor exclusions, genotype/truth integrity checks and explicit failure reporting; successful workflow execution alone is not proof of valid scientific results.
- Leave already-running experiments and their frozen producers/protocols unchanged. Apply this preference to the next version or new experiment rather than interrupting a run to retrofit orchestration.
- Small standalone checks and unit tests may run directly; they do not require a Nextflow wrapper.
- After an analysis completes and its published results pass independent verification, remove that run's disposable Nextflow work/cache data so disk usage does not accumulate. Scope cleanup to the exact run in its dedicated launch directory; retain results, models, input/source hashes, execution reports and failure evidence. Keep failed/incomplete runs resumable until explicitly recovered or retired. Never use global cache pruning or delete shared scientific inputs.

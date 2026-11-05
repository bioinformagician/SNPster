#!/bin/sh
exec /opt/snpster-nextflow-metrics/lib/ld-linux-x86-64.so.2 \
    --library-path /opt/snpster-nextflow-metrics/lib \
    /opt/snpster-nextflow-metrics/bin/ps "$@"

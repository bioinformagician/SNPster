#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: bash prepare_nextflow_metrics.sh ABSOLUTE_OUTPUT" >&2; exit 1; fi
output="$(realpath -m "$1")"
project="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -d "$output" ]]; then
    test -s "$output/assets.sha256"
    sha256sum --check "$output/assets.sha256"
    exit
fi
mkdir -p "$output/bin" "$output/lib"
image='sha256:d2ce7a8fad9afdd38c943c7bb3d0495d718a32305e1aac96e22d9f32e42a769b'
helper="$(docker create --network none "$image" /bin/true)"
trap 'docker rm "$helper" >/dev/null' EXIT
docker cp -L "$helper:/usr/bin/ps" "$output/bin/ps"
for library in libprocps.so.8 libc.so.6 libsystemd.so.0 liblzma.so.5 libzstd.so.1 \
    liblz4.so.1 libcap.so.2 libgcrypt.so.20 libgpg-error.so.0; do
    docker cp -L "$helper:/lib/x86_64-linux-gnu/$library" "$output/lib/$library"
done
docker cp -L "$helper:/lib64/ld-linux-x86-64.so.2" "$output/lib/ld-linux-x86-64.so.2"
cp "$project/pedigree_benchmark/scripts/nextflow_ps.sh" "$output/ps"
chmod 755 "$output/ps" "$output/bin/ps" "$output/lib/ld-linux-x86-64.so.2"
printf '%s\n' "$image" > "$output/source_image.txt"
sha256sum "$output/ps" "$output/bin/ps" "$output/lib/"* "$output/source_image.txt" > "$output/assets.sha256"
sha256sum --check "$output/assets.sha256"

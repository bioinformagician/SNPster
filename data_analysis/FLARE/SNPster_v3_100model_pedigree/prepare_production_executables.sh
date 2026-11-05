#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then echo "Usage: bash prepare_production_executables.sh ABSOLUTE_OUTPUT" >&2; exit 1; fi
output="$(realpath -m "$1")"
if [[ -d "$output" ]]; then
    test -x "$output/plink"
    test -x "$output/plink2"
    sha256sum --check "$output/assets.sha256"
    exit
fi
mkdir -p "$output"
image='sha256:cd7907c83ad8d41017ce1e235172ace6c5c00d04565a11f75c0b185260c7a74f'
helper="$(docker create --network none --entrypoint /usr/bin/sha256sum "$image" \
    /app/dependencies/plink /app/dependencies/plink2)"
trap 'docker rm "$helper" >/dev/null' EXIT
for binary in plink plink2; do
    docker cp -L "$helper:/app/dependencies/$binary" "$output/$binary"
    chmod 755 "$output/$binary"
done
docker start -a "$helper" > "$output/original.sha256"
if [[ "$(docker inspect --format '{{.State.ExitCode}}' "$helper")" != 0 ]]; then
    echo "Could not obtain executable hashes from the original pinned image." >&2; exit 1
fi
while read -r expected filename; do
    actual="$(sha256sum "$output/${filename##*/}")"
    if [[ "${actual%% *}" != "$expected" ]]; then
        echo "Extracted executable differs from original pinned image: $filename" >&2; exit 1
    fi
done < "$output/original.sha256"
printf '%s\n' "$image" > "$output/source_image.txt"
sha256sum "$output/plink" "$output/plink2" "$output/original.sha256" \
    "$output/source_image.txt" > "$output/assets.sha256"
sha256sum --check "$output/assets.sha256"

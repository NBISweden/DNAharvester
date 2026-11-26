#!/usr/bin/env bash
# Download an oras:// image into the Nextflow Apptainer cache, skipping Apptainer's architecture
# check. Needed for the arm64 G-ALN image when Nextflow runs on an x86 node: `apptainer pull`
# refuses it there, but Nextflow uses the image if it is already in its cache.
#
# Usage: fetch_oras_image.sh oras://<registry>/<repo>:<tag> <cache_dir>
set -euo pipefail

ref=${1#oras://}                                           # registry/repo:tag
registry=${ref%%/*}; repo=${ref#*/}; tag=${repo##*:}; repo=${repo%:*}
out=${2:?usage: $0 oras://<registry>/<repo>:<tag> <cache_dir>}/$(echo "$ref" | tr '/:' '--').img
url=https://$registry/v2/$repo

# Anonymous pull token, then the digest of the image (the single layer of the manifest)
realm=$(curl -sI "$url/manifests/$tag" | sed -nE 's/.*realm="([^"]+)".*/\1/p')
token=$(curl -s "$realm?service=$registry&scope=repository:$repo:pull" | python3 -c 'import sys,json; print(json.load(sys.stdin)["token"])')
digest=$(curl -s -H "Authorization: Bearer $token" -H 'Accept: application/vnd.oci.image.manifest.v1+json' \
    "$url/manifests/$tag" | python3 -c 'import sys,json; print(json.load(sys.stdin)["layers"][0]["digest"])')

curl -sL -H "Authorization: Bearer $token" -o "$out" "$url/blobs/$digest"
echo "${digest#sha256:}  $out" | sha256sum -c --quiet || { rm -f "$out"; exit 1; }
echo "Saved: $out"

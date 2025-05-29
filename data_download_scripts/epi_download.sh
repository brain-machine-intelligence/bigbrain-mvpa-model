#!/usr/bin/env bash
set -e
echo "[EPI] Downloading multi-voxel pattern data from OSF..."
curl -L "https://files.osf.io/v1/resources/2gyue/providers/osfstorage/6836f7ed8ff21163f3d6a5c3/?zip=" -o epi.zip   
mkdir -p data
unzip -q epi.zip -d data/
rm epi.zip
echo "[EPI] Done. Files are in data/subj_masked_EPI/"

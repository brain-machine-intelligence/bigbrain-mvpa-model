#!/usr/bin/env bash
set -e
echo "[ZIP] Downloading compressed demo data..."
curl -L "https://osf.io/v27cy/download" -o data.zip
mkdir -p data
unzip -q data.zip -d data/
rm data.zip
echo "[ZIP] Done."

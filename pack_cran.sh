#!/bin/bash
# packs the current git HEAD (or any other supplied git head) into a tarball
# suitable for CRAN submission. Name comes from git describe as the latest tag.

N=EmbedSOM
HEAD="${1:-HEAD}"
VERTAG=$(git describe --tags --no-abbrev "${HEAD}")
VER=${VERTAG#v}
ARCHIVE=${N}_${VER}.tar.gz

git archive --format=tar --prefix="${N}/" "${HEAD}" | \
tar -f - \
	--delete "${N}/pack_cran.sh" \
	--delete "${N}/vignettes/.gitignore" \
	--delete "${N}/src/.clang-format" | \
gzip -c > "${ARCHIVE}"

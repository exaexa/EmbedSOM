#!/bin/bash
# packs the current git HEAD (or any other supplied git head) into a tarball
# suitable for CRAN submission. Name comes from git describe as the latest tag.

N=EmbedSOM
HEAD="${1:-HEAD}"
VERTAG=$(git describe --tags --no-abbrev "${HEAD}")
VER=${VERTAG#v}
ARCHIVE=${N}_${VER}.tar.gz

TMPDIR=".tmpbuild-$$"

mkdir ${TMPDIR} || exit 1

git archive --format=tar --prefix="${N}/" "${HEAD}" \
| tar f - \
	--delete "${N}/pack_cran.sh" \
	--delete "${N}/vignettes/.gitignore" \
	--delete "${N}/src/.clang-format" \
	--delete "${N}/src/style.sh" \
| tar xf - -C ${TMPDIR}

R CMD build ./${TMPDIR}/${N}/ --compact-vignettes

rm -fr ${TMPDIR}

R CMD check ${ARCHIVE}

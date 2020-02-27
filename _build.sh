#!/bin/sh

set -ev

Rscript -e "pkgdown::build_site()"

set -o errexit -o nounset

[ -z "${GH_TOKEN}" ] && exit 0
[ "${TRAVIS_BRANCH}" != "master" ] && exit 0

git config --global user.email "psolymos@gmail.com"
git config --global user.name "Peter Solymos"

git clone -b gh-pages https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git stats-output
cd stats-output
cp -r ../_stats/* ./
git add --all *
git commit -m "Update stats (${TRAVIS_BUILD_NUMBER})" || true
git push -q origin gh-pages

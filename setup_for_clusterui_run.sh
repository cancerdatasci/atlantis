set -xe
mkdir ./tRpkg
export R_LIBS=./tRpkg:/data2/Rpackages
R CMD INSTALL . -l ./tRpkg

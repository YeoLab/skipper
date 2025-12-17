#!/bin/bash

#  Install UMICollapse v1.0.0
echo "Installing UMICollapse"

curl -L -O https://github.com/Daniel-Liu-c0deb0t/UMICollapse/archive/refs/tags/v1.0.0.tar.gz
tar -xzf v1.0.0.tar.gz
rm *gz

cd UMICollapse-1.0.0

mkdir lib
cd lib

curl -L -O https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar
curl -L -O https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar

echo "Setup complete."

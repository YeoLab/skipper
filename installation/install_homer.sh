#!/bin/bash

# Install HOMER
echo "Installing HOMER"

CURRENT=$(pwd)

# Install HOMER
cd $CURRENT
mkdir homer && cd ./homer \
    && curl -O http://homer.ucsd.edu/homer/configureHomer.pl \
    && perl configureHomer.pl -install \
    && cd $CURRENT

HOMER_PATH="$CURRENT/homer/bin"

echo "export PATH=$HOMER_PATH:\$PATH" >> ~/.bashrc

echo "Setup complete."

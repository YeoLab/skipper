#!/bin/bash

#  Install samtools
echo "Installing samtools"

CURRENT=$(pwd)

# Install samtoolswhich 
cd $CURRENT
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
tar -xjf samtools-1.17.tar.bz2 && \
cd samtools-1.17 && \
./configure --prefix=$CURRENT/samtools-1.17-build --without-curses --disable-bz2 && \
make CPPFLAGS=-I$HOME/opt/bzip2/include LDFLAGS="-L$HOME/opt/bzip2/lib -Wl,-R$HOME/opt/bzip2/lib" && \
make install && \
cd $CURRENT && \
rm -rf $CURRENT/samtools-1.17 $CURRENT/samtools-1.17.tar.bz2

SAMTOOLS_PATH="$CURRENT/samtools-1.17-build/bin"

echo "export PATH=$SAMTOOLS_PATH:\$PATH" >> ~/.bashrc

echo "Setup complete."

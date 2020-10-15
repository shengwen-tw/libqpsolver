#!/bin/bash

#download and run mkl installer
wget https://registrationcenter-download.intel.com/akdlm/irc_nas/tec/16903/l_mkl_2020.3.279.tgz
tar -zxvf l_mkl_2020.3.279.tgz
cd l_mkl_2020.3.279
sudo ./install_GUI.sh

#append LD_LBRARY_PATH export command in .bashrc
INSTALL_PATH=/opt/intel/mkl
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'$INSTALL_PATH'/lib/intel64/' >> ~/.bashrc

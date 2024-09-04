#!/bin/bash

#if [[ $# -eq 0 ]] ; then
#    echo "Initialize the virtual environment for gwskysim"
#    exit 0
#fi

#create Virtual Environment
current_dir=`pwd`
parent_dir="$(dirname "$current_dir")"

venv_dir=${parent_dir}/venv-gwskysim-py3
/usr/bin/python3.8 -m venv --without-pip ${venv_dir}
echo $venv_dir
source ${venv_dir}/bin/activate
curl https://bootstrap.pypa.io/get-pip.py | python
deactivate
source ${venv_dir}/bin/activate

pip install --upgrade pip setuptools wheel
pip install -r ${current_dir}/requirements.txt
deactivate

#Create a setup script
echo "source ${venv_dir}/bin/activate" > ${current_dir}/setup-gwskysim.sh
echo "export PATH=${current_dir}/bin:\$PATH" >> ${current_dir}/setup-gwskysim.sh
echo "export PYTHONPATH=${current_dir}:\$PYTHONPATH" >> ${current_dir}/setup-gwskysim.sh
echo "echo">> ${current_dir}/setup-gwskysim.sh
echo "echo \"gwskysim setup done!\"">> ${current_dir}/setup-gwskysim.sh
echo "echo">> ${current_dir}/setup-gwskysim.sh


echo "Initialization done!. To setup the package:"
echo "source setup-gwskysim.sh"

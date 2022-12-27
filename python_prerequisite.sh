if command -v wget >/dev/null 2>&1; then
  wget -O ./requirements.txt https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/requirements.txt
else
  if command -v curl >/dev/null 2>&1; then
    curl -o ./requirements.txt https://raw.githubusercontent.com/ShuXingYu94/NGS_Process_Sample/master/requirements.txt
  else
    echo "Wget and Curl are not installed. Please install."
    exit 0
  fi
fi

if command -v pip3 >/dev/null 2>&1; then
  pip3 install -r requirements.txt
else
  if command -v pip >/dev/null 2>&1; then
    pip install -r requirements.txt
  else
    echo "Something wrong with pip/pip3. Please install python packages with other software like conda."
    exit 0
  fi
fi

rm ./requirements.txt
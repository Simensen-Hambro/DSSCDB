#!/bin/bash
printenv > /etc/environment

cd /code/
mkdir -p media/spreadsheets
mkdir -p media/molecules
cd /code/src/
sleep 7
python manage.py migrate
if [ "$DEBUG" == "False" ]
then
  python manage.py collectstatic --noinput
fi

source activate rdkitenv
gunicorn DSSCDB.wsgi:application --log-file "-" --access-logfile "-" --log-level INFO -w 2 -b :8000

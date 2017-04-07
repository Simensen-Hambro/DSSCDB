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

echo "##### STARTING UWSGI CRON AS DAEMON #####"
uwsgi -d /dev/null --module=DSSCDB.cron_wsgi:application --master --pidfile=/tmp/project-cron.pid --http=0.0.0.0:8001 --processes=5 --harakiri=20 --max-requests=5000 --vacuum

echo "##### UWSGI MAIN #####"
uwsgi --module=DSSCDB.wsgi:application --master --pidfile=/tmp/project-master.pid --http=0.0.0.0:8000 --processes=5 --harakiri=20 --max-requests=5000 --vacuum

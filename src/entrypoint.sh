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
uwsgi -d /dev/stderr --module=DSSCDB.cron_wsgi:application --pidfile=/tmp/project-cron.pid --master --logto /dev/stdout --http=0.0.0.0:8001 --processes=1 --harakiri=300 --enable-threads --vacuum

echo "##### UWSGI MAIN #####"
uwsgi --module=DSSCDB.wsgi:application --master --pidfile=/tmp/project-master.pid --logto /dev/stderr --http=0.0.0.0:8000 --processes=5 --harakiri=300 --max-requests=5000 --vacuum
#/opt/conda/envs/rdkitenv/bin/python -u /opt/.pycharm_helpers/pydev/pydevd.py --multiprocess --qt-support --port 53300 --file /code/src/manage.py runserver 8000

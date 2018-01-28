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

python manage.py loaddata email_templates.json

echo "##### STARTING UWSGI CRON AS DAEMON #####"
uwsgi -d /dev/null --module=DSSCDB.cron_wsgi:application --pidfile=/tmp/project-cron.pid --master --http=0.0.0.0:8001 --processes=1 --harakiri=300 --enable-threads --vacuum

echo "##### UWSGI MAIN #####"
uwsgi --module=DSSCDB.wsgi:application --master --pidfile=/tmp/project-master.pid --http=0.0.0.0:8000 --processes=5 --harakiri=300 --max-requests=5000 --vacuum

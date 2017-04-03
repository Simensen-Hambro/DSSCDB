#!/bin/bash
printenv > /etc/environment
cd /code/
sleep 7
python manage.py migrate
if [ "$DEBUG" == "False" ]
then
  python manage.py collectstatic --noinput
fi

exec uwsgi --chdir=./ --module=DSSCDB.wsgi:application --master --pidfile=/tmp/project-master.pid \
--http=0.0.0.0:8000 --processes=5 --max-requests=5000 --vacuum

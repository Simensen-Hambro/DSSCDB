rm db.sqlite3
find . -iname "migrations" -exec rm -rf {} \;
python manage.py makemigrations DSSCDB dye usermanagement
python manage.py migrate
python manage.py loaddata test_data.json


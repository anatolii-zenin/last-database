(was stored in DESY gitlab for 2 years)
# LAST Polarisation Database

A tool for storing all observation data in a single database with two tables:
1. "diffphot" for raw diffphot output data
2. "poldata" for processed polarisation data

## Getting started

Requires mysql-server package. <br>
Mysql must have user "root" with password "admin" and a database "mydatabase". <br>
To install mysql and to change the password one has to run a series of commands: <br>
```
$sudo apt-get install mysql-server
$sudo service mysql stop
$sudo killall mysqld
$sudo chown -R mysql:mysql /var/lib/mysql
$sudo mysqld_safe --skip-grant-tables #(might have to use chown for another folder if you get an error)
$sudo mysql -u root -p <br>
```
in mysql: <br>
```
mysql> update mysql.user set authentication_string=password('NEWPASSWORD') where user='root'; #(replace "NEWPASSWORD" with your password)
mysql> flush privileges;
mysql> quit <br>
```
back in shell: <br>
```
$sudo killall mysqld
$sudo service mysql stop
$sudo service mysql start
$mysql -u root -p #(your new password should work now)
$mysql > CREATE DATABASE mydatabase;
```
## Current functionality:

1. Read a folder with .instr files into a dataframe and save it in the database
2. Calculate polarisation data from the database and save it in a new dataframe and save it in the database
3. Retrieve diffphot data from a database
4. Retrieve a polarisation dataframe from a database

## How to use:

### Reading .instr files and filling a diffphot database:
db = DatabaseHandler("./data/db_test/", "diffphot")

###

## Issues:

1. Diffphot files do not contain any information about the coordinates of the source 
2. The databases can not be appended automatically
3. Database edits should edit input files?
4. Mysql configuration can be problematic

## To do:

1. Add an "append database" method with overlap checks?
2. Add a method to remove data from the database (and from the original files?)


import helpers
import openroad as ord
import odb
import sys

from ctypes import CDLL
libc = CDLL('libc.so.6')

db = ord.get_db()

odb.read_lef(db, "data/gscl45nm.lef")
odb.read_def(db, "data/design.def")

db_file = "results/export.db"

write_result = odb.write_db(db, db_file)

if write_result != 1:
    print("FAIL: Write DB failed")
    exit(1)

new_db = odb.dbDatabase_create()
odb.read_db(new_db, db_file)

# db_diff returns True is there is a difference.
if odb.db_diff(db, new_db):  
    print("FAIL: Differences found between exported and imported db")
    exit(1)

# If we do not flush the buffer for C/C++ then "pass" will not be printed
# last and the test will fail
libc.fflush(None)
print("pass")
exit(0)

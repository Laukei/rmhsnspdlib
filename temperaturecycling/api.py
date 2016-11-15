#
# A library for working around temperature-cycling of a cryostat with a non-indefinite hold time
#
# By Rob Heath rob@robheath.me.uk 15/11/2016
#

import os
import sqlite3
import datetime

folder = os.path.split(os.path.abspath(__file__))[0]

db_filename = os.path.join(folder,'record.db') #filename for database
db_schema = os.path.join(folder,'schema.sql') #schema for db creation
get_old = 60 #time in s after which values should be considered old and invalid

class RecordKeeper:
	def __init__(self,**kwargs):
		self.db_filename = db_filename if 'db_filename' not in kwargs else kwargs['db_filename']
		self.db_schema = db_schema if 'db_schema' not in kwargs else kwargs['db_schema']
		self.get_old = get_old if 'get_old' not in kwargs else kwargs['get_old']
		self.old_override = False if 'old_override' not in kwargs else kwargs['old_override']
		self._check_db()

	def _check_db(self):
		db_is_new = not os.path.exists(self.db_filename)
		with sqlite3.connect(self.db_filename) as conn:
			if db_is_new:
				with open(self.db_schema, 'rt') as f:
					schema = f.read()
				conn.executescript(schema)

	def add_temperature(self,temperature):
		self._check_db()
		with sqlite3.connect(self.db_filename) as conn:
			conn.execute('INSERT INTO temperature (base_temp) VALUES (?)',(temperature,))

	def get_temperature(self):
		self._check_db()
		with sqlite3.connect(self.db_filename) as conn:
			cur = conn.cursor()
			cur.execute('SELECT id, added_time, base_temp FROM temperature ORDER BY id DESC LIMIT 1')
			result = cur.fetchone()
		return result

	def can_continue(self,threshold_temp):
		self.t = self.get_temperature()
		if self.t != None:
			self.t_since = ((datetime.datetime.now()-datetime.datetime.strptime(self.t[1],'%Y-%m-%d %H:%M:%S')).total_seconds())
			self.t_temp = self.t[2]
			if ((self.t_since <= self.get_old) and self.t_temp <= threshold_temp) or self.old_override == True:
				return True
		elif self.t == None and self.old_override == True:
			return True
		return False

if __name__ == '__main__':
	rk = RecordKeeper()
	rk.add_temperature(3.5)
	print(rk.can_continue(3.4))

from datetime import datetime
from pathlib import Path
import logging
import sys

class Tee:
	def __init__(self, *files):
		self.files = files

	def write(self, obj):
		for f in self.files:
			f.write(obj)

	def flush(self):
		for f in self.files:
			f.flush()

def print_log(log_file_name):
	current_time = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")
	log_file_name = f"{log_file_name}"
	logging.basicConfig(filename=log_file_name, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
	log_file = open(log_file_name, 'a')
	sys.stdout = Tee(sys.stdout, log_file)
	print("-"*100)
	print("-"*100)
	print(f">>> Date and Time: {current_time}")
	print(f">>> logfile name: {log_file_name}")
	print("-"*100)
	print("-"*100)
	return

if __name__ == '__main__':
	logs_folder = Path("logs")
	logs_folder.mkdir(parents=True, exist_ok=True)
	log_file_name = f"./{logs_folder}/print.log"
	print_log(log_file_name)

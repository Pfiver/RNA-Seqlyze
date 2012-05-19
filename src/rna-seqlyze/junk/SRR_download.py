#!/usr/bin/python
# encoding: utf-8

import ftplib, re, urllib

host = 'ftp-trace.ncbi.nlm.nih.gov'

#basedir = '/sra/sra-instant/reads/ByStudy/sra/SRP/SRP001/SRP001461/'

ftp = ftplib.FTP(host)
ftp.login('anonymous', 'patrick.pfeifer@students.fhnw.ch')
ftp.cwd(basedir)

def dirwalk(ftp, callback):    
	print 'Path:', ftp.pwd()
  
	dirs = []
	def cb(line):
		if line.startswith("d"):
			dir = re.split(r' *', line, maxsplit=8)[8]
			if dir not in ('.', '..'):
				dirs.append(dir)
	ftp.dir("", cb)

	try:
		for item in dirs:
			ftp.cwd(item)
			print 'Changed to', ftp.pwd()
			callback(ftp)
			dirwalk(ftp, callback)
			ftp.cwd('..')
	except Exception, e:
		print item, e

def collect_urls(ftp):
	files = []
	def cb(line):
		if line.startswith("-"):
			files.append(re.split(r' *', line, maxsplit=8)[8])
	ftp.dir("", cb)

	wd = ftp.pwd()
	for fil in files:
		urls.append("ftp://" + host + wd + "/" + fil)

urls = []

dirwalk(ftp, collect_urls)

for url in urls:
	print "Fetching", url
	def cb(blocks, bloksize, total):
		print "Status: %d of %d\r" % (blocks*blocksize, total)
	filename, message = urllib.urlretrieve(url, url.split("/")[-1], cb)
	print
	print "Done: %s" % filename
	print message.headers
	print message.fp.read()

"""pyCasJobs module.

This module provides a class for accessing the SDSS CasJobs server in
a programmatic manner.

Now discovered the method used (pretending to be a real user, storing
cookies) is unnecessary as CasJobs provides an API.
E.g., see http://casjobs.sdss.org/CasJobs/services/jobs.asmx?op=UploadData

Example usage:
>>> import pyCasJobs
>>> cas = pyCasJobs.CasJobs('username', 'password')
>>> cas.import_table('mytable.csv', 'mytable')

The current functionality is limited only to importing tables.
However, it should be straightforward to implement the additional elements of CasJobs,
such as SQL queries, in a similar manner.

Requires a non-standard library module: poster, see http://atlee.ca/software/poster/.

Created on Jul 15, 2010

@author: Steven Bamford
"""

import urllib2, urllib, cookielib
from poster.encode import multipart_encode
import poster.streaminghttp

class CasJobs:
    """Programmatic access to the SDSS CasJobs server.
    
    A class which enables programmatic access to the Sloan Digital Sky Survey
    CasJobs server, by wrapping the usual browser-based web interface.

    Keyword arguments to constructor:
    username -- CasJobs username
    password -- CasJobs password

    """
    def __init__(self, username, password):
        self.username = username
        self.password = password
        self.loginurl = 'http://casjobs.sdss.org/CasJobs/login.aspx'
        self.importurl = 'http://casjobs.sdss.org/CasJobs/TableImport.aspx'
        opener = poster.streaminghttp.register_openers()
        opener.add_handler(urllib2.HTTPCookieProcessor(cookielib.CookieJar()))
        self.login()

    def login(self):
        params = urllib.urlencode({'userid': self.username, 'password': self.password})
        response = urllib2.urlopen(self.loginurl, params)
        self.page = response.read()
        response.close()

    def import_table(self, filename, tablename, tableexists=False, format='text',
                     ntries=3):
        """Upload a local file into CasJobs MyDB.
    
        Note that CasJobs has rather stringent limits to the size of file
        which can be uploaded.

        Keyword arguments:
        filename -- filename of the table to upload
        tablename -- name of the table to create/append to in CasJobs MyDB
        tableexists -- if 'False' create a new table, if 'True' append to existing table
        format -- 'text': space/comma/tab separated text file
                  'votable': XML VOTable
                  'dataset': MS DataSet XML

        """
        if tableexists:
            tableTypeDDL = 1
        else:
            tableTypeDDL = 0
        if format == 'text':
            DataType = 0
        elif format == 'votable':
            DataType = 1
        elif format == 'dataset':
            DataType = 2
        datagen, headers = multipart_encode({'tableTypeDDL': tableTypeDDL,
                                             'NameBox': tablename, 'tnameDDL': tablename,
                                             'DataType': DataType, 'importType': 1,
                                             'DataFormat': 3, 'sigmaBox': 5, 'DataBox': '',
                                             'httpBox': open(filename, 'rb')})
        request = urllib2.Request(self.importurl, datagen, headers)
        tries = 0
        while True:
            try:
                tries += 1
                response = urllib2.urlopen(request)
                break
            except urllib2.URLError:
                if tries > ntries:
                    raise
        self.page = response.read()
        response.close()

if __name__ == '__main__':
    pass

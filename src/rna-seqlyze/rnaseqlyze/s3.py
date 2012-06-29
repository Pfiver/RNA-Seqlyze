"""
Upload files to Amazon S* using the s3cmd tools

The access credentials must be configured in ~/.s3cfg .
The s3cmd creates this file interactively with the --configure option.
"""

#: The s3 bucket name in "S3Uri" format
base_uri = "s3://biocalc/"

import os

# import this to fix a circular import dependency problem in s3cmd ...
import S3.Exceptions

from S3.S3 import S3
from S3.Config import Config
from S3.S3Uri import S3Uri
from S3.SortedDict import SortedDict

# same "interface" like rnaseqlyze.galaxy
def upload(fileobj, filename):

    cfg = Config()

    cfg.read_config_file(os.path.join(os.getenv("HOME"), ".s3cfg"))
    cfg.progress_meter = False
    cfg.acl_public = True

    s3 = S3(cfg)

    headers = SortedDict(ignore_case = True)
    headers["x-amz-acl"] = "public-read"
    headers["x-amz-storage-class"] = "REDUCED_REDUNDANCY"

    remote_uri = S3Uri(base_uri + filename)

    fileobj.seek(0,2) # seek to end
    size = fileobj.tell()
    fileobj.seek(0) # seek to start

    response = s3.send_file_multipart(fileobj, headers, remote_uri, size)

    assert response['status'] == 200

    return remote_uri.public_url()

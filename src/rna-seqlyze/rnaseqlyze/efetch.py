from Bio import Entrez

import rnaseqlyze

nc_db = "nuccore"
gb_type = "gb"
gb_mode = "text"

Entrez.email = rnaseqlyze.entrez_email

def get_nc_id(accession):
    handle = Entrez.esearch(db=nc_db, term=accession + "[Accession]")
    id_list = Entrez.read(handle)["IdList"]
    if len(id_list) != 1:
        raise Exception("unexpected reply from Entrez: id_list: %s" % id_list)
    return id_list[0]

def fetch_nc_gb(gb_id, out_file):
    handle = Entrez.efetch(db=nc_db, id=gb_id, rettype=gb_type, retmode=gb_mode)
    from shutil import copyfileobj as copy
    copy(handle, out_file)

protocol = 'http'
hostname = 'archaea.ucsc.edu'
path_reset = '/cgi-bin/cartReset'
path_tracks = '/cgi-bin/hgTracks'
custom_text = 'hgt.customText'

custom_track_url = "%(protocol)s://%(hostname)s%(path_tracks)s" \
                        "?db=%%(db)s&%(custom_text)s=%%(track_url)s" % globals()

# http://archaea.ucsc.edu/cgi-bin/hgTracks?db=sulSol1&hgt.customText=

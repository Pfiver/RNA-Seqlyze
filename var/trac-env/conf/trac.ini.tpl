# -*- coding: utf-8 -*-

[tractab]
names = BuildBot, TestCoverage
urls = /@@WWWBASE_@@buildbot/console, /@@WWWBASE_@@trac/chrome/site/files/coverage/index.html
perms = WIKI_VIEW, WIKI_VIEW

[navadd]
add_items = app, apidocs
app.title = Application
app.url = /@@WWWBASE_@@rna-seqlyze
app.target = mainnav
apidocs.title = Apidocs
apidocs.url = /@@WWWBASE_@@trac/chrome/site/files/apidoc/index.html
apidocs.target = mainnav

[attachment]
max_size = -1
max_zip_size = -1
render_unsafe_content = false

[browser]
color_scale = True
downloadable_paths = /trunk, /branches/*, /tags/*
hide_properties = svk:merge
intermediate_color = 
intermediate_point = 
newest_color = (255, 136, 136)
oldest_color = (136, 136, 255)
oneliner_properties = trac:summary
render_unsafe_content = false
wiki_properties = trac:description

[changeset]
max_diff_bytes = 10000000
max_diff_files = 0
wiki_format_messages = true

[components]
tracopt.ticket.deleter = enabled
footnotemacro.macro.footnotemacro = enabled
tracopt.ticket.commit_updater.committicketreferencemacro = enabled
tracopt.ticket.commit_updater.committicketupdater = enabled
tracopt.versioncontrol.git.* = enabled

[header_logo]
alt = Header Logo
height = -1
link = /
src = site/img/RNA-Seqlyze.png
width = -1

[inherit]
htdocs_dir = 
plugins_dir = 
templates_dir = 

[logging]
log_file = @@WORKDIR@@-dev/trac.log
log_level = INFO
log_type = file

[mainnav]
browser.label = Source Code
browser.href = /browser/src?rev=master
search = disabled
tickets.href = /report/2

[metanav]
about = disabled

[milestone]
stats_provider = DefaultTicketGroupStatsProvider

[mimeviewer]
max_preview_size = 262144
mime_map = text/x-trac-wiki:trac, text/x-dylan:dylan, text/x-idl:ice, text/x-ada:ads:adb
tab_width = 8
treat_as_binary = application/octet-stream, application/pdf, application/postscript, application/msword,application/rtf,

[notification]
admit_domains = 
always_notify_owner = false
always_notify_reporter = false
always_notify_updater = true
ambiguous_char_width = single
batch_subject_template = $prefix Batch modify: $tickets_descr
email_sender = SmtpEmailSender
ignore_domains = 
mime_encoding = none
sendmail_path = sendmail
smtp_always_bcc = 
smtp_always_cc = 
smtp_default_domain = 
smtp_enabled = false
smtp_from = trac@localhost
smtp_from_author = false
smtp_from_name = 
smtp_password = 
smtp_port = 25
smtp_replyto = trac@localhost
smtp_server = localhost
smtp_subject_prefix = __default__
smtp_user = 
ticket_subject_template = $prefix #$ticket.id: $summary
use_public_cc = false
use_short_addr = false
use_tls = false

[project]
admin = 
admin_trac_url = .
descr = 
footer = Visit the Trac open source project at<br /><a href="http://trac.edgewall.org/">http://trac.edgewall.org/</a>
icon = site/img/RNA-Seqlyze-icon.ico
name = RNA-seqlyze
url = 

[query]
default_anonymous_query = status!=closed&cc~=$USER
default_query = status!=closed&owner=$USER
items_per_page = 100
ticketlink_query = ?status=!closed

[report]
items_per_page = 100
items_per_page_rss = 0

[revisionlog]
default_log_limit = 100
graph_colors = ['#cc0', '#0c0', '#0cc', '#00c', '#c0c', '#c00']

[roadmap]
stats_provider = DefaultTicketGroupStatsProvider

[search]
min_query_length = 3

[svn]
branches = trunk, branches/*
tags = tags/*

[ticket]
commit_ticket_update_notify = true
commit_ticket_update_envelope =
commit_ticket_update_check_perms = false
commit_ticket_update_commands.refs = <ALL>
commit_ticket_update_commands.close = close closed closes fix fixed fixes
default_cc = 
default_component = Software
default_description = 
default_keywords = 
default_milestone = Version 1.0
default_owner = patrick
default_priority = medium
default_resolution = fixed
default_severity = 
default_summary = 
default_type = defect
default_version = 1.0
max_comment_size = 262144
max_description_size = 262144
preserve_newlines = default
restrict_owner = false
workflow = ConfigurableTicketWorkflow

[ticket-workflow]
accept = new,assigned,accepted,reopened -> accepted
accept.operations = set_owner_to_self
accept.permissions = TICKET_MODIFY
leave = * -> *
leave.default = 1
leave.operations = leave_status
reassign = new,assigned,accepted,reopened -> assigned
reassign.operations = set_owner
reassign.permissions = TICKET_MODIFY
reopen = closed -> reopened
reopen.operations = del_resolution
reopen.permissions = TICKET_CREATE
resolve = new,assigned,accepted,reopened -> closed
resolve.operations = set_resolution
resolve.permissions = TICKET_MODIFY

[timeline]
abbreviated_messages = True
changeset_collapse_events = false
changeset_long_messages = false
changeset_show_files = 0
default_daysback = 30
max_daysback = 90
newticket_formatter = oneliner
ticket_show_details = false

[trac]
auth_cookie_lifetime = 0
auth_cookie_path = 
authz_file = 
authz_module_name = 
auto_preview_timeout = 2.0
auto_reload = False
backup_dir = db
base_url = http://@@HOSTNAME@@/@@WWWBASE_@@trac/
check_auth_ip = false
database = @@TRAC_DB@@
debug_sql = False
default_charset = utf-8
default_date_format = iso8601
default_dateinfo_format = relative
default_language = 
default_timezone = GMT +2:00
genshi_cache_size = 128
htdocs_location = 
ignore_auth_case = false
jquery_location = 
jquery_ui_location = 
jquery_ui_theme_location = 
mainnav = app, apidocs, BuildBot, TestCoverage, wiki, roadmap, tickets, newticket, browser, timeline
metanav = login, logout, prefs, help
mysqldump_path = mysqldump
never_obfuscate_mailto = false
permission_policies = DefaultPermissionPolicy, LegacyAttachmentPolicy
permission_store = DefaultPermissionStore
pg_dump_path = pg_dump
repository_sync_per_request =
resizable_textareas = true
secure_cookies = False
show_email_addresses = false
show_ip_addresses = false
timeout = 20
use_base_url_for_redirect = False

[repositories]
.type = git
.dir = @@TOPDIR@@/.git
.url =

[versioncontrol]
allowed_repository_dir_prefixes = 

[wiki]
ignore_missing_pages = false
max_size = 262144
render_unsafe_content = false
safe_schemes = cvs, file, ftp, git, irc, http, https, news, sftp, smb, ssh, svn, svn+ssh
split_page_names = false


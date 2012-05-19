# gbMembers.sql was originally generated by the autoSql program, which also 
# generated gbMembers.c and gbMembers.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

# This sql was hacked to insert the datetime object, autoSql could not
# do that and added auto_increment to the idx field and index with
# userName

#UCSC Genome Browser members
DROP TABLE IF EXISTS gbMembers;
CREATE TABLE gbMembers (
    idx int unsigned NOT NULL auto_increment,	# auto-increment unique ID
    userName varchar(255) NOT NULL,	# Name used to login
    realName varchar(255) NOT NULL,	# Full name
    password varchar(255) NOT NULL,	# Encrypted password
    email varchar(255) NOT NULL,	# Email address
    lastUse DATETIME NOT NULL, # Last date the user log in/log out/change password
    newPassword varchar(255) NOT NULL, # Password generated for the mail-a-new-password feature
    newPasswordExpire DATETIME NOT NULL, # Expiration date of the new password generated
    dateActivated DATETIME NOT NULL, # Date the account activated via email
    emailToken varchar(255) NOT NULL, # Security token used in the email to the user
    emailTokenExpires DATETIME NOT NULL,	# Expiration date of the emailToken
    passwordChangeRequired char(1) NOT NULL default 'N',	# Password change required?
    accountActivated char(1) NOT NULL default 'N',	# Account activated? Y or N
              #Indices
    PRIMARY KEY(idx),
    INDEX(userName)
);
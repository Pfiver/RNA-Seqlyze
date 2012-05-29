# history

histfile=~/.pbh
if [ -d $(dirname $histfile) -a "$HISTFILE" != $histfile ]
then
	HISTSIZE=1000000
        HISTFILE=$histfile
        HISTIGNORE="&:sync"
        HISTTIMEFORMAT="%Y-%m-%d %H:%M:%S  "
        shopt -s histappend

        if [ -f "$HISTFILE" ] && [ $(wc -l < $HISTFILE) -ge $HISTSIZE ]
        then
                echo "warning: HISTFILE ('$HISTFILE') reached size limit of $HISTSIZE lines"
        fi
fi

# other stuff

set -PC

CDPATH=:~

xpath() {
        [[ ${!1} =~ (^|:)$2(:|$) ]] || eval ${1}=\$2:\$$1
}
xpath PATH ~/.local/bin

alias vi=vim

alias cp='cp -i'
alias mv='mv -i'
alias rm='rm -i'
alias grep='grep --color=auto'
alias ls='ls --color=auto --show-control-chars --group-directories-first'

alias gs='git status'

eval "$(lesspipe)"
eval "$(dircolors -b)"

stty -ixon -ixoff

export VISUAL=vim EDITOR=vim # PYTHONSTARTUP=~/.config/python/startup.py

pushd() {
	if [ $# = 0 ]
	then
		if [ ${#DIRSTACK[@]} -le 2 ]
		then
			command pushd "$PWD" > /dev/null
		else
			cd ~-0
			command pushd -n -0 > /dev/null
		fi
	elif [ -d "$1" ]
	then
		command pushd -n "$(readlink -f "${1%/}")" > /dev/null
		cd ~1
	else
		echo no such directory: "'$1'" >&2
	fi
}
popd() {
	command popd "$@" > /dev/null
}

PROMPT_COMMAND='[ ${#DIRSTACK[@]} -gt 1 ] && { echo -e "\033[01;31m --\033[0m"; dirs -v | tail -n+2; }'


PS1='\[\033[;92m\]\u@\h\[\033[00m\]:\[\033[;94m\]\w\[\033[00m\]\$ '

alias egit="git --git-dir=/home/pfeifer/.git --work-tree=/home/pfeifer"

bt() {
	bt=/home/pfeifer/data/bt
	cd $bt
	. bash-env
}

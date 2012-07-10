# tmux
######

tmux="$HOME/.local/bin/tmux -S /tmp/tmux-pfeifer"

if ! [ $TMUX ]
then
	$tmux display 2> /dev/null && exec $tmux attach
	TERM=xterm-256color exec $tmux
fi

# history
#########

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
                echo "warning: HISTFILE ('$HISTFILE') \
			reached size limit of $HISTSIZE lines"
        fi
fi

# other stuff
#############

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

alias gd='git diff'
alias gs='git status'

alias egit="git --git-dir=/home/pfeifer/.git --work-tree=/home/pfeifer"

eval "$(lesspipe)"
eval "$(dircolors -b)"

stty -ixon -ixoff

export VISUAL=vim EDITOR=vim PYTHONSTARTUP=~/.pythonstartup.py

PS1='\[\033[;92m\]\u@\h\[\033[00m\]:\[\033[;94m\]\w\[\033[00m\]\$ '
PROMPT_COMMAND='echo -ne "\033]0;${USER}@$(hostname -s): ${PWD/$HOME/~}\007"'

find-py()
{
	find "$@" \( -name junk -o -name setup.py \) -prune -o \
		\( -name "*.pt" -o \( -name "*.py" \( ! -name __init__.py \
		-o -size +1 -exec sed -n '1{/^__/{q1}}' {} \; \) \) \) -print
}

# viall
#######
# edit all python, javascript, page-templte and lessccss files
# in the RNA-Seqlyze source tree
viall ()
{
    cd $rnas_topdir/src

    not_starts_with__="-exec sed -n 1{/^__/{q1}} {} ;"

    vi $(find rna-seqlyze* \(                           \
            -name junk -o -name setup.py                \
        \) -prune -o \(                                 \
            -name "rnaseqlyze*.js"                      \
            -o -name "rnaseqlyze*.less"                 \
            -o -name "*.pt"                             \
            -o -name "*.py" \(                          \
                -not -name __init__.py                  \
                -o -size +1 $not_starts_with__ \)       \
        \) -print |
                        tr '/-' '[]' | sort | tr '[]' '/-'
    ) 
}


# fix pushd
###########
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
RED=$'\033[01;31m'
BOLD=$'\033[01;1m'
NORM=$'\033[0m'

PROMPT_COMMAND+='; [ "0${#DIRSTACK[@]}" -gt 1 ] && \
			{ echo -e " $RED--$NORM"; dirs -v | tail -n+2; }'

# rna-seqlyze
#############

bt() {
	. ~/data/rna-seqlyze/bash-env

	dirs -c
	for ((i=${#rnas_pkgdirs[@]}-1; i>0; i--))
	do
		pkg=${rnas_pkgdirs[i]##*-}
		pushd "${rnas_pkgdirs[i]}/rnaseqlyze/$pkg"
	done
	pushd "${rnas_pkgdirs[0]}/rnaseqlyze/core"
	pushd "${rnas_pkgdirs[0]}/rnaseqlyze"

	bt=$rnas_topdir
	wd=$rnas_workdir

	cd "$rnas_topdir"
}

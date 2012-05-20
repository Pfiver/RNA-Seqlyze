# history

histfile=~/.pbh
if [ -d $(dirname $histfile) -a $HISTFILE != $histfile ]
then
        HISTSIZE=1000000
        HISTFILE=$histfile
        HISTIGNORE="&:sync"
        HISTFILESIZE=1000000
        HISTTIMEFORMAT="%Y-%m-%d %H:%M:%S  "
        shopt -s histappend

        if [ -f "$HISTFILE" ] && [ $(wc -l < $HISTFILE) -ge $HISTFILESIZE ]
        then
                echo "warning: HISTFILE ('$HISTFILE') reached size limit of $HISTFILESIZE lines"
        fi

        history -c
        history -r
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

eval "$(lesspipe)"
eval "$(dircolors -b)"

stty -ixon -ixoff

export VISUAL=vim EDITOR=vim # PYTHONSTARTUP=~/.config/python/startup.py

PS1='\[\033[;92m\]\u@\h\[\033[00m\]:\[\033[;94m\]\w\[\033[00m\]\$ '

bt() { bt=/home/pfeifer/data/bt; cd $bt; . bash-env; }

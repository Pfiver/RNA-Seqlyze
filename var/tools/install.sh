#!/bin/bash

(cd inotify-tools-3.14 &&
	./configure --prefix=$HOME/.local && make all install)

(cd libevent-2.0.19-stable
	&& ./configure --prefix=$HOME/.local && make all install)

(cd tmux-1.6 &&
	LIBEVENT_CFLAGS= LIBEVENT_LIBS=-levent \
	PKG_CONFIG=true LDFLAGS=-L$HOME/.local/lib \
	CFLAGS=$(echo -I$HOME/.local/include{,/ncurses}) \
	./configure --prefix=$HOME/.local &&
		LD_RUN_PATH=$HOME/.local/lib make all install)

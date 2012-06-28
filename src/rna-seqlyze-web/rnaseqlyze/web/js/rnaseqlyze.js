/*
 * RNA-seqlyze javascript routines
 */

$(function() {

    // global variables and helpers
    // ----------------------------

    // ___ el ___
    //  from http://joestelmach.github.com/laconic/
    window.el = $.el;

    // log.info() and log.debug()
    window.log = {
        'info': function() {
            console.log.apply(console, arguments);
        },
        'debug': function() {}
    }
    if (rnaseqlyze_debug)
        window.log.debug = window.log.info;

    // page initializaion
    // ------------------

    // use bootstrap's
    // "scrollspy" plugin
    //  -- patched version - see https://github.com/twitter/bootstrap/pull/3829
    $(window).scrollspy({
//            offset: 200,
            wrap: $('#wrap')[0],
    });

    $('#content .page-header').sticky({
            topSpacing: 0,
    });

});

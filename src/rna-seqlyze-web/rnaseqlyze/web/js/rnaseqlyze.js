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


    // Generally useful stuff
    // ----------------------

    // based on http://stackoverflow.com/a/4673436

    String.prototype.format = function() {
        var i = 0; args = arguments;
        return this.replace(/{}/g, function() {
            return args[i++];
        });
    };

});

// http://stackoverflow.com/a/7531350
jQuery.fn.extend({
    scrollTo: function () {
        var offset = jQuery(this).offset();
        if (!offset) return;
        var x = offset.top - 100;
        jQuery('html,body').animate({scrollTop: x}, 100);
    },
});

/*
 * RNA-seqlyze 'analysis' view javascript
 */

// backbone.js Models
// ------------------
//
// -> http://backbonejs.org/#Model

// this analysis
window.Analysis = Backbone.Model.extend({
    urlRoot: "../rest/analyses",
    initialize: function () {
        this.files = new WorkdirListing();
        this.files.analysis = this;

        // "cascade": update files list
        this.bind("change", function (self) {
            // FIXME: maybe update this less often
            //        in production, e.g. every second time,
            self.files.fetch();
        });
    }
});
// a model for the files
window.WorkdirFile = Backbone.Model.extend({
    defaults: {
        path: null,
    }
});
// and for a collection of files
window.WorkdirListing = Backbone.Collection.extend({
    model: WorkdirFile,
    url: function () {
        return this.analysis.url() + "/files";
    },
});

// note1:
//
//  Concerning the above code:
//  It might have been simpler to work the files list right into the
//  analysis model on the server and stick to one model here.
//  But then again, there is no harm in doing it like this, because
//  now the files list is more independent and could for
//  example also be displayed on a page of its own.


/* note2:
 *
 *  In the code below,
 *
 *  "el"           is defined in rnaseqlyze.js as "$.el", which is defined
 *                 in laconic.js - see  http://joestelmach.github.com/laconic/
 *
 *  "this.$el"
 *  "this.el"      are the view's (jQuery wrapped) DOM element in the
 *                 backbone.js architecture - see http://backbonejs.org/#View-el
 *
 *  "render().el"  is also the view's "el" and works because we always
 *                 "return this;" from render() - see http://backbonejs.org/#View-render
 */

// Two Views showing different details about the analysis
// ------------------------------------------------------

// These render the "Processing" and "Results" section on the
// analysis page. The Processing view is displayed above the Results view.
//
// -> http://backbonejs.org/#View

// The "Processing" section
window.ProcessingView = Backbone.View.extend({

    initialize: function () {
        this.model.bind("change", this.change, this);
    },
    change: function (model, value, options) {

        // just re-render the whole thing for now
        this.$el.empty();
        this.render();

        // remove the busy indicator when finished
        if (model.get('finished'))
            $('#spinner').remove();

        // make scrollspy refresh it's coordinates
        // because the page size has very likely changed
        // DON'T FIXME: this sadly gets called twice now
        //              since we have two views...
        //              no harm in it! :-)
        $(window).scrollspy('refresh');
    },
    render: function () {
        // toJSON doesn't really do much besides turning
        // the model.attributes into a useable object
        // see http://backbonejs.org/#Model-toJSON
        var analysis = this.model.toJSON();

        this.$el.append(el.div(
            el.h3("Input analysis")
        ,
            analysis.inputfile_name ?
                el.p("Type of input: ",
                     analysis.inputfile_type ?
                        el.strong(analysis.inputfile_type) :
                        el.span("not detected"))
                :
                null
        ,
            analysis.inputfile_header ?
                el.p("First lines in input data: ",
                     el.pre(analysis.inputfile_header))
                :
                null

        )).append(el.div(
            el.h3("Working Directory")
        ,
            el.div(new WorkdirView({model: this.model.files}).render().el)
        ));

        return this;
    },
});

// A View displaying the list of files associated with
// this analysis available on the server (log files, mostly).
// This is currently rendered inside the "Processing" section above.
window.WorkdirView = Backbone.View.extend({
    initialize: function () {
        this.model.bind("reset", this.reset, this);
    },
    reset: function (model, value, options) {
        this.$el.empty();
        this.render();
    },
    render: function () {
        this.$el.append(el.ul.apply(el,
            _(this.model.models).map(function (model) {
                return new WorkdirFileView({model: model}).render().el;
            })
        ));
        return this;
    },
});
// An View, that renders one file
window.WorkdirFileView = Backbone.View.extend({
    el: "<li>",
    render: function (model) {
        var file = this.model.toJSON();
        var href = _id + '/files/' + file.path;
        this.$el.html(el.a({href: href}, file.path));
        return this;
    }
});

// The "Results" section
window.ResultsView = Backbone.View.extend({
    initialize: function () {
        this.model.bind("change", this.change, this);
    },
    change: function (model, value, options) {
        this.$el.empty();
        this.render();
        $(window).scrollspy('refresh');
    },
    render: function () {
        var analysis = this.model.toJSON();
        if (analysis.hg_url)
            this.$el.append(
                el.ul(
                    el.li(
                        el.a({href: analysis.hg_url},
                             "Link to BAM Track in UCSC Browser"))));
        return this;
    },
});

// Initialization
// --------------

$(document).ready(function () {

    // the id of the displayed analysis
    _id = _(window.location.pathname.split('/')).last();

    // create a backbone.js Model
    // with an associated Collection
    analysis = new Analysis({
        id: _id,
    });

    // create two backbone.js views for the
    // analysis, render and insert them into the DOM
    $('#processing-area').html(
        new ProcessingView({model: analysis}).render().el
    );
    $('#results-area').html(
        new ResultsView({model: analysis}).render().el
    );

    // see what's going on with backbone.js
    if (rnaseqlyze_debug) {
        analysis.bind("all", function (event) {
               log.debug("analysis", event);
        });
        analysis.files.bind("all", function (event) {
                log.debug("analysis.files", event);
        });
    }

    // update the models until the analysis is finished
    var update = function () {
        // check at the beginning and not at the end
        // because the fetch() calls are asynchronous
        if (analysis.attributes.finished)
            return;
        analysis.fetch();
        log.debug(analysis);
        window.setTimeout(update, 7000); // re-update in 7 seconds
    }
    update();
});

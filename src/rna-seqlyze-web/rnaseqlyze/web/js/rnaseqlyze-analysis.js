/*
 * RNA-seqlyze 'analysis' view javascript
 */

// backbone.js Models
// ------------------
//
// -> http://backbonejs.org/#Model

// The Analysis
window.Analysis = Backbone.Model.extend({
    urlRoot: "../rest/analyses",
    initialize: function () {
        this.files = new DataDirListing();
        this.files.analysis = this;

        this.stage_logs = new StageLogList();
        this.stage_logs.analysis = this;

        // "cascade": update files list
        this.bind("change:data_dir_state", function (self) {
            self.files.fetch({add: true});
        });
        // "cascade": update stage_logs
        this.bind("change:stage_logs_state", function (self) {
            var len = self.stage_logs.size();
            self.stage_logs.fetch({
                add: true,
                // The last stage_log is the current one and updates frequently.
                // It's id stays the same though and which causes backbone.js
                // to regard it as a duplicate and drop it. But a copy of the
                // ajax response is passed to the success callback. So we pick
                // the changed log text from there and fire a "change" event
                // by set()ting the 'text' attribute of the affected model.
                success: function (stage_logs, rsp) {
                    if (len)
                        stage_logs.at(len-1).set('text', rsp[len-1].text);
                },
            });
        });
    }
});
// log output of one stage
window.StageLog = Backbone.Model.extend({
    defaults: {
        stage: null,
        text: null,
    },
    idAttribute: "id",
});
// log output of all stages
window.StageLogList = Backbone.Collection.extend({
    model: StageLog,
    url: function () {
        return this.analysis.url() + "/logs";
    },
});
// a model for the files
window.DataDirFile = Backbone.Model.extend({
    defaults: {
        path: null,
    },
    idAttribute: "path",
});
// and for a collection of files
window.DataDirListing = Backbone.Collection.extend({
    model: DataDirFile,
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
        this.stage_logs = (
            new StageLogListView({model: this.model.stage_logs}).render().el);
    },
    change: function (model, value, options) {

        // just re-render the whole thing for now
        this.$el.empty();
        this.render();

        // remove the busy indicator when finished
        if (model.get('finished'))
            $('#spinner').remove();

        // make scrollspy refresh it's coordinates
        // because the page size has likely changed
        $(window).scrollspy('refresh');
    },
    render: function () {
        // toJSON doesn't really do much besides turning
        // the model.attributes into a useable object
        // see http://backbonejs.org/#Model-toJSON
        var analysis = this.model.toJSON();

        this.$el.append(
            el.h2("Processing")
        );

        this.$el.append(el.div(
            el.h3("Input Check")
        ,
            analysis.inputfile_uploaded ?
                el.p("Type of input: ",
                     analysis.inputfile_type ?
                        el.strong(analysis.inputfile_type) :
                        el.span("not detected"))
                :
                null
        ,
            analysis.inputfile_header ?
                el.p("First read in input data: ",
                     el.pre(analysis.inputfile_header))
                :
                null

        ));

        this.$el.append(el.div(
            el.h3("Stage Logs")
        ,
            this.stage_logs
        ));

        if (analysis.error)
            this.$el.append(
                el.div({class: "alert alert-error"},
                       el.h4({class: "alert-heading"},
                             "An error occured while analyzing the data"),
                       analysis.error));

        return this;
    },
});

// The monospaced stage log blocks
window.StageLogListView = Backbone.View.extend({
    initialize: function () {
        this.model.bind("add", this.add, this);
        if (!this.model.analysis.get('finished')) {
            this.model.analysis.bind("change:finished",
                                     this.analysis_change, this);
        }
    },
    analysis_change: function (model, value, options) {
        if (value) {
            this.$el.contents().find("pre")
                    .last().css('background-color', '');
        }
    },
    add: function (model) {
        this.$el.append(
            new StageLogView({model: model}).render().el
        );
        if (!this.model.analysis.get('finished')) {
            this.$el.contents().find("pre")
                .not(':last').css('background-color', '')
                .prevObject.last().css('background-color', '#ddf');
            this.$el.contents().last().scrollToBottom();
        }
        $(window).scrollspy('refresh');
    },
});

// _One_ monospaced stage log block
window.StageLogView = Backbone.View.extend({
    initialize: function () {
        this.model.bind("change", this.change, this);
    },
    change: function (model, options) {
        this.$el.children("pre").text(model.get('text'));
        this.$el.scrollToBottom();
        $(window).scrollspy('refresh');
    },
    render: function () {
        var log = this.model.toJSON();
        this.$el.attr('id', log.stage);
        this.$el.append(el.h4(log.stage));
        this.$el.append(el.pre(log.text));
        return this;
    },
});

// The "Results" section
window.ResultsView = Backbone.View.extend({
    initialize: function () {
        this.model.bind("change", this.change, this);
        this.model.files.bind("add", this.fileadd, this);
    },
    change: function (model, value, options) {
        this.$el.empty();
        this.render();
        $(window).scrollspy('refresh');
    },
    fileadd: function () {
        var augmented_gb = this.model.files.find(function (file) {
            return file.get('path').match(/augmented\.gb$/);
        });
        if (!this.augmented_gb) {
            this.augmented_gb = augmented_gb;
            this.$el.empty();
            this.render();
        }
    },
    render: function () {
        var analysis = this.model.toJSON();
        if (analysis.hg_url || this.augmented_gb) {
            this.$el.append(el.h2("Results"));
            var $ul = $(el.ul())
            this.$el.append($ul[0]);
            if (this.augmented_gb) {
                var href = _id + '/files/' + this.augmented_gb.get('path');
                $ul.append(el.li(
                    el.a({href: href},
                         "Augmented Genbank File")));
            }
            if (analysis.hg_url) {
                $ul.append(el.li(
                    el.a({href: analysis.hg_url},
                         "Link to custom tracks in UCSC browser"),
                    el.p("It might take a minute until the tracks become " +
                         "available.", el.br(),
                         "As soon as the last few items ",
                         el.a({href: galaxy_history_url}, "here"),
                         " turn green it should work.")));
            }
        }
        return this;
    },
});

// A View displaying the list of files associated with
// this analysis available on the server (log files, mostly).
// This is currently rendered inside the "Processing" section above.
window.DataDirView = Backbone.View.extend({
    initialize: function () {
        this.model.bind("reset", this.reset, this);
        this.model.bind("add", this.add, this);
    },
    reset: function (model, value, options) {
        this.$el.empty();
        this.render();
    },
    render: function () {
        this.$el.append(el.h2("Data Directory"));
        var ul = el.ul();
        this.$ul = $(ul);
        this.$el.append(ul);
        $(window).scrollspy('refresh');
        return this;
    },
    add: function(model) {
        this.$ul.append(
            new DataDirFileView({model: model}).render().el
        );
        $(window).scrollspy('refresh');
    },
});
// An View, that renders one file
window.DataDirFileView = Backbone.View.extend({
    el: "<li>",
    render: function (model) {
        var file = this.model.toJSON();
        var href = _id + '/files/' + file.path;
        this.$el.html(el.a({href: href}, file.path));
        return this;
    }
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
    $('#processing').html(
        new ProcessingView({model: analysis}).render().el
    );
    $('#results').html(
        new ResultsView({model: analysis}).render().el
    );
    $('#datadir').html(
        new DataDirView({model: analysis.files}).render().el
    );

/*
    // uncomment this to see what's going on in backbone.js
    if (rnaseqlyze_debug) {
        analysis.bind("all", function (event) {
            log.debug("analysis", arguments);
        });
        analysis.files.bind("all", function (event) {
            log.debug("analysis.files", arguments);
        });
        analysis.stage_logs.bind("all", function (event) {
            log.debug("analysis.stage_logs", arguments);
        });
    }
 */

    // update the models until the analysis is finished
    var update = function () {
        // check at the beginning and not at the end
        // because the fetch() calls are asynchronous
        if (analysis.attributes.finished)
            return;
        analysis.fetch();
//        log.debug("analysis.fetch()");
        window.setTimeout(update, 7000); // re-update in 7 seconds
    }
    update();
});

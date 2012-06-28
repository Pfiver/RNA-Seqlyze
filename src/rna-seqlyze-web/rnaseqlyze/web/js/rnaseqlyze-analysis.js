/*
 * RNA-seqlyze 'analysis' view javascript
 */

// Custom Configuration
_.templateSettings = { interpolate : /\{(.+?)\}/g };

// Globals
_id = _(window.location.pathname.split('/')).last();
console = {'log': function(){}}

// Models
window.Analysis = Backbone.Model.extend({
    urlRoot: "../rest/analyses",
});
window.File = Backbone.Model.extend({
    idAttribute: "path",
});
window.FileList = Backbone.Collection.extend({
    model: File,
    url: "../rest/analyses/" + _id + "/files",
});

// Views
window.AnalysisView = Backbone.View.extend({
 
    initialize: function () {
        this.model.bind("change", this.change, this);
        this.model.bind("reset", this.render, this);
        this.model.bind("all", function(event) { console.log(event); });
    },

    change: function(model, value, options) {

        if (!this.model.attributes.inputfile_header)
            this.model.attributes.inputfile_header = '';

        console.log(this.model);
        console.log(this.model.toJSON());
        console.log("json:" + this.model.toJSON());
        console.log("mvo: " + model + ", " + value + ", " + options);

         $(this.el).html(
             _.template($('#tpl-input-type').html())(this.model.toJSON()));

         $(this.el).html(
             _.template($('#tpl-results-list-hg-url').html())(this.model.toJSON()));
    },
    render: function () {
        console.log(this.model);
//        $(this.el).html(this.template(this.model.toJSON()))
        return this;
    },
});    

window.FileListView = Backbone.View.extend({
 
    initialize: function () {
        this.model.bind("change", this.render, this);
        this.model.bind("reset", this.render, this);
        this.model.bind("all", function(event) { console.log(event); });
    },
 
    render: function () {
        console.log(this.model);
        _(this.model.models).each(function (file) {
            $(this.el).append(new FileListItemView({model: file}).render().el);
        }, this);
        return this;
    },
});    

window.FileListItemView = Backbone.View.extend({
 
    template: _.template($('#tpl-logfile-list-item').html()),

    initialize: function () {
        this.model.bind("change", this.render, this);
        this.model.bind("add", this.render, this);

        this.model.bind("all", function(event) { console.log(event); });
    },

    render: function () {
        console.log(this.model);
        $(this.el).html(this.template(this.model.toJSON()));
        return this;
    },
});

// Initialization
$(document).ready(function() {

    analysis = new Analysis({id: _id});
    file_list = new FileList();

    analysis_view = new AnalysisView({model: analysis});
    file_list_view = new FileListView({model: file_list});

    analysis.fetch();
    file_list.fetch();

    console.log(analysis);
    console.log(file_list);

    $('#input-type').html(analysis_view.render().el);
    $('#logfile-list').html(file_list_view.render().el);

/*
    var interval = window.setInterval(function() {
        fileList.fetch();
        return true;
    }, 7000);
 */

});

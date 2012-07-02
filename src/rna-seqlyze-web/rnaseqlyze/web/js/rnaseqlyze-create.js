/*
 * RNA-seqlyze javascript routines
 *
 *     for the "create" page
 */

$(document).ready(function() {

    /*
     * Toggle input type
     */

    $('#input_type_radio').click(function(event) {
        if ($(event.target).hasClass("srr")) {
            // user chose the "Data File" option
            $('#sra-controls').hide();
            $('#srr-controls').show();
        } else if ($(event.target).hasClass("sra")) {
            // user chose the "SRR Identifier" option
            $('#srr-controls').hide();
            $('#sra-controls').show();
        } // else
            // what the...
    });
    $('#input_type_radio .srr').click();

    /*
     * discretionary #pairendlenControls
     */

    function maybe_show_pairendlen_controls() {
        if ($('#pairendedInput').attr('checked'))
            $('#pairendlenControls').show();
        else
            $('#pairendlenControls').hide();
    }

    maybe_show_pairendlen_controls();
    $('#pairendedInput').change(maybe_show_pairendlen_controls);


    /*
     * Toggle organism input type
     */

    $('#org_type_radio').click(function(event) {
        if ($(event.target).hasClass("title")) {
            // user chose "Title"
            $('#genbankfile-controls').hide();
            $('#org_title-controls').show();
        } else if ($(event.target).hasClass("file")) {
            // user chose "Genbank File"
            $('#org_title-controls').hide();
            $('#genbankfile-controls').show();
        } // else
            // what the ...
    });
    $('#org_type_radio .title').click();

    /*
     * Organism input autocompletion
     */

    var organisms = new Array();

    $.ajax({
        url: "rest/organisms",
        dataType: "json",
        success: function(data) {
            // No idea yet what to do with those that have multiple
            // accessions listed in 'genome' -- filter them out here for now
            //
            // see ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
            //
            _(data).each(function(org) {
                if (org.acc.indexOf(",") >= 0)
                    return;
                organisms.push("{} ({}/{})".format(
                                org.title, org.db, org.acc)); }); },
    });

    $('#organismInput').typeahead({
        source: organisms,
    });

    /*
     * plupload -- from src/plupload/examples/custom.html
     */

    var uploads = {
        'inputfile': {},
        'genbankfile': {}
    };

    var context = function(name) {

        var self = this;

        var options = {
            url:      'upload',
            runtimes:      'html5,gears,flash,silverlight,browserplus,html4',
            browse_button:    name + '_browse',
            drop_element:       name + '_progress',
            multipart_params:     { 'type': name,
                                      'session': upload_session },
            flash_swf_url:          path_js + '/plupload.flash.swf',
            silverlight_xap_url:      path_js + '/plupload.silverlight.xap',
        };

        var events = {
            // 'Init': function(up, params) {
            //     $('#' + name + '_progress .filestatus').text(
            //                     "Current runtime: " + params.runtime);
            // },
            'FilesAdded': function(up, up_files) {
                // remove all other files already present
                // plupload features multiple files in one widget
                // we have two widgets and one name per widget
                up.splice();
                self.active = true;
                log.debug("FilesAdded", up.files, up_files);
                $('#' + name + '_progress .filestatus').text(
                    up_files[0].name +
                        ' (' + plupload.formatSize(up_files[0].size) + ')');
            },
            'UploadComplete': function(up, up_files) {
                self.complete = true;
                for (nam in uploads)
                    if (uploads[nam].active)
                        if(!uploads[nam].complete)
                            return;
                log.debug("UploadComplete", "go");
                $('#create_form').submit();
            },
            'UploadProgress': function(up, up_file) {
                $('#' + name + '_progress .bar').css(
                            "width", up_file.percent + '%');
                // $('#' + name + '_progress .filestatus').text(
                //                             up_file.percent + '%');
            },
        };

        this.active = false;
        this.complete = false;

        var up = this.up = new plupload.Uploader(options);

        $('#' + name + '_progress').click(function() {
            $('#' + name + '_browse').click();
        });

        $('#create_form_submit').click(function() {
            for (nam in uploads)
                if (uploads[nam].active)
                    { up.start(); return false; }
        });

        for (x in events)
            up.bind(x, events[x]);

        log.debug("debug");

        up.init();
    };

    for (name in uploads)
        uploads[name] = new context(name);
});

// vim: et:sw=4

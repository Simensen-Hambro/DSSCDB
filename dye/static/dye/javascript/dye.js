$('.pagination-button').click(function (e) {
    e.preventDefault();
    var $form = $("#search-form");
    $form.attr("action", $(this).attr("href"));
    $form.submit();
});

$('.file-upload-form').submit(function (e) {
    $('.loader').css("visibility", "visible");
    $('.loader-text').css("display", "block");
});

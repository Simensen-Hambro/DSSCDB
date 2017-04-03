$(".clickable-row").click(function () {
    window.document.location = $(this).data("href");
});

$(".clickable-cell").click(function () {
   window.document.location = $(this).parent().data("href");
});

function PlotlyHoverGroups(){};

PlotlyHoverGroups.HoverGroups = {};

PlotlyHoverGroups.resetHoveredPoints = function(group_id) {
    const plots = PlotlyHoverGroups.HoverGroups[group_id];
    plots.forEach(plot => {
        Plotly.Fx.unhover(plot);
    });
};

PlotlyHoverGroups.findPointIndex = function(plot_id, point_id) {
    const plot = document.getElementById(plot_id);
    const data = plot.data;
    const index = data[0].customdata.findIndex(x => x === point_id);
    return index;
};

PlotlyHoverGroups.setHoveredPoint = function(group_id, point_id, calling_plot_id = null) {
    let plot_ids = PlotlyHoverGroups.HoverGroups[group_id];
    if (calling_plot_id !== null) {
        plot_ids = plot_ids.filter(plot_id => plot_id !== calling_plot_id);
    }
    plot_ids.forEach(plot_id => {
        const point_index = PlotlyHoverGroups.findPointIndex(plot_id, point_id);
        Plotly.Fx.hover(plot_id, [{ curveNumber: 0, pointNumber: point_index }]);
    });
};

Shiny.addCustomMessageHandler('add-plotly-hover-group', function(message) {
    console.log(message);
    PlotlyHoverGroups.HoverGroups[message.id] = message.members;
});


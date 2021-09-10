figure(1)
clf
hold on
c = ["#5e81ac", "#8fbcbb", "#bf616a", "#d08770", "#ebcb8b", "#a3be8c", "#b48ead"];
cid  = [2 1 3 3 4 4 5 5];
show = [1 1 1 0 1 0 1 0];
name = ["Cio", "Bio", "Dio", "", "Eio", "", "Fio", ""];
for i = [2, 1, 3:8]
    x = wpddatasets(:, 2*i-1);
    y = wpddatasets(:, 2*i);
    idx = find(x < 0);
    if ~isempty(idx)
        x(idx) = x(idx) + 180;
        x = [x(1:idx(end)); NaN; x(idx(end)+1:end)];
        y = [y(1:idx(end)); NaN; y(idx(end)+1:end)];
    end
    idx = find(x > 180);
    if ~isempty(idx)
        x(idx) = x(idx) - 180;
        x = [x(1:idx(1)-1); NaN; x(idx(1):end)];
        y = [y(1:idx(1)-1); NaN; y(idx(1):end)];
    end
    if show(i); bool = 'on'; else; bool = 'off'; end
    plot(x, y, "linewidth", 2, "color", c(cid(i)), 'HandleVisibility', bool, "DisplayName", name(i))
end
grid on
xlabel('Solar Phase Angle (deg)')
ylabel('Time of Flight(days)')
legend('FontSize', 12, 'Location', 'best')
hold off
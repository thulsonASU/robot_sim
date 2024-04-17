% function event to stop integration and final time
function [value, isterminal, direction] = events(t, x, tf)
    value = tf - t;     % When value is 0, an event is triggered
    isterminal = 1;     % Stop the integration
    direction = 0;      % The zero can be approached from either direction
end
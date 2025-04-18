function myPrint(app, formatStr, varargin)

   msg = sprintf(formatStr, varargin{:});
    lines = splitlines(msg);        
    lines = cellstr(lines);        
    app.ConsoleOutput.Value = [app.ConsoleOutput.Value; lines];
    scroll(app.ConsoleOutput, 'bottom')
    drawnow;
end

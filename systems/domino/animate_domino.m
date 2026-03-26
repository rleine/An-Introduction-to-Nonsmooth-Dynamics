function animate_domino(q, n, a, b, saveVideo, skip)

    if nargin < 5, saveVideo = false; end
    if nargin < 6, skip = 0; end

    Nt = size(q, 2);        % number of time steps

    figure('Color','w');
    hold on;
    axis equal;
    axis off;               
    box off;

    % Axis limits 
    xmin = min(q(1:3:end, :), [], 'all') - 2*b;
    xmax = max(q(1:3:end, :), [], 'all') + 2*b;
    ymax = max(q(2:3:end, :), [], 'all') + 3*b;

    axis([xmin xmax -0.05 ymax]);
    
     % Draw the ground
    groundHeight = b/2;
    groundX = [xmin, xmax, xmax, xmin];
    groundY = [0, 0, -groundHeight, -groundHeight];

    patch(groundX, groundY, [0.7 0.7 0.7], ...
      'EdgeColor','none');

    set(gcf,'Position',[100 100 600 125])  


    % Predefine rectangle in body frame K (centered at origin)
    rect_K = [
        -a, -b;
         a, -b;
         a,  b;
        -a,  b;
        -a, -b
    ]';

    % Create patch objects for dominoes
    h = gobjects(n,1);
    for j = 1:n
        h(j) = patch('XData', [], 'YData', [], ...
                     'FaceColor', [0.8 0.3 0.3], ...
                     'EdgeColor', 'k');
    end

    if saveVideo 
        filename = fullfile(pwd, 'domino_simulation.avi');
        v = VideoWriter(filename);
        v.FrameRate = 50;
        v.Quality = 70;
        open(v);
    end

    x = q(1:3:end,:);
    y = q(2:3:end,:);
    phi = q(3:3:end,:);
    c = cos(phi); s = sin(phi);

    % Animation loop
    for i = 1:skip:Nt

        for j = 1:n

            % Rotation matrix
            A_IK = [c(j,i) -s(j,i);
                    s(j,i)  c(j,i)];

            % Rotate and translate rectangle
            rect_I = A_IK * rect_K + kron(ones(1,5),[x(j,i);y(j,i)]);

            % Update patch
            set(h(j), 'XData', rect_I(1,:), ...
                      'YData', rect_I(2,:));
        end

        drawnow;

        if saveVideo 
            frame = getframe(gcf);
            writeVideo(v, frame);
        end

    end

    if saveVideo
        close(v);
    end
end

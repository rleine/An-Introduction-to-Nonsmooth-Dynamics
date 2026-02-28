function animate_superball(q, R, saveVideo, skip, pause_time)

    if nargin < 3, saveVideo = false; end
    if nargin < 4, skip = 1; end
    if nargin < 5, pause_time = 0; end

    Nt = size(q,2);

    x   = q(1,:);
    y   = q(2,:);
    phi = q(3,:);

    c = cos(phi);
    s = sin(phi);

    %% Figure setup
    figure('Color','w');
    hold on
    axis equal
    axis off
    box off

    % Drawing the floor
    floormin=-0.5;
    floorwidth=2;
    floorheight=floorwidth;
    heightmin=-1;
    axis([floormin floormin+floorwidth heightmin heightmin+floorheight]);
    plot([floormin,floormin+floorwidth],[0 0],'k','LineWidth',2); 


    % Drawing the ceiling
    ceilingmin= 0.3;
    ceilingwidth=1;
    l = 1;
    plot([ceilingmin,ceilingmin+ceilingwidth],l+[0 0],'k','LineWidth',2);

    %% Precompute geometry in body frame

    th = linspace(0,2*pi,100);
    circle_K = [R*cos(th); R*sin(th)];

    % Cross lines in body frame
    line1_K = [ R  0;
               -R  0]';

    line2_K = [ 0  R;
                0 -R]';

    %% Initial transformation
    A_IK = [c(1) -s(1);
             s(1)  c(1)];

    circle_I = A_IK*circle_K + [x(1); y(1)];
    line1_I  = A_IK*line1_K  + [x(1); y(1)];
    line2_I  = A_IK*line2_K  + [x(1); y(1)];

    %% Create graphics objects
    h_circle = plot(circle_I(1,:), circle_I(2,:), 'b','LineWidth',1.5);

    h_line1  = plot(line1_I(1,:), line1_I(2,:),   'r','LineWidth',1.5);

    h_line2  = plot(line2_I(1,:), line2_I(2,:),   'r','LineWidth',1.5);

    %% Video writer
    if saveVideo
        filename = fullfile(pwd,'superball_animation.avi');
        v = VideoWriter(filename);
        v.FrameRate = 50;
        v.Quality   = 80;
        open(v);
    end

    %% Animation loop
    for i = 1:skip:Nt

        % Rotation matrix
        A_IK = [c(i) -s(i);
                s(i)  c(i)];

        % Transform geometry
        circle_I = A_IK*circle_K + [x(i); y(i)];
        line1_I  = A_IK*line1_K  + [x(i); y(i)];
        line2_I  = A_IK*line2_K  + [x(i); y(i)];

        % Update graphics
        set(h_circle,'XData',circle_I(1,:), ...
                     'YData',circle_I(2,:));

        set(h_line1,'XData',line1_I(1,:), ...
                    'YData',line1_I(2,:));

        set(h_line2,'XData',line2_I(1,:), ...
                    'YData',line2_I(2,:));

        drawnow;
        pause(pause_time);

        if saveVideo
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end

    if saveVideo
        close(v);
    end

end

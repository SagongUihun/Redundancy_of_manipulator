clear all;
close all;
%% Make pose, velocity, accleration of End-Effector
EE_start=[0.4, 0.0, 0.5];
EE_goal=[0.4, 0, 0.4];

tf = 7;
ta = tf/3;

vel_max = (EE_goal - EE_start) / (tf-ta);
accel_max = vel_max / ta;

i=0;
time_step=0.1;
for t=0:time_step:tf
    i=i+1;

    if t < ta
        accel_EE(i,:) = accel_max;
        vel_EE(i,:) = accel_max*t;
        pos_EE(i,:) = EE_start + accel_max/2*t*t;

    elseif t < tf-ta
        accel_EE(i,:) = 0;
        vel_EE(i,:) = vel_max;
        pos_EE(i,:) = EE_start + vel_max/2*ta + vel_max*(t-ta);

    else
        accel_EE(i,:) = -accel_max;
        vel_EE(i,:) = -accel_max*(t-tf);
        pos_EE(i,:) = EE_goal - accel_max/2*(t-tf)*(t-tf);
    end
end

%% Visualization of EE
% for j=1:size(pos_EE,1)
%     figure(1)
%     subplot(3,1,1)
%     plot(pos_EE(:,1));
%     grid on
%     ylabel('position');
% 
%     subplot(3,1,2)
%     plot(vel_EE(:,1));
%     grid on
%     ylabel('velocity');
% 
%     subplot(3,1,3)
%     plot(accel_EE(:,1));
%     grid on
%     xlabel('time(sec)');
%     ylabel('acceleration');
% 
% end
% 
% for j=1:size(pos_EE,1)
%     figure(2)
%     subplot(3,1,1)
%     plot(pos_EE(:,2));
%     grid on
%     ylabel('position');
% 
%     subplot(3,1,2)
%     plot(vel_EE(:,2));
%     grid on
%     ylabel('velocity');
% 
%     subplot(3,1,3)
%     plot(accel_EE(:,2));
%     grid on
%     xlabel('time(sec)');
%     ylabel('acceleration');
% 
% end
% 
% for j=1:size(pos_EE,1)
%     figure(3)
%     subplot(3,1,1)
%     plot(pos_EE(:,3));
%     grid on
%     ylabel('position');
% 
%     subplot(3,1,2)
%     plot(vel_EE(:,3));
%     grid on
%     ylabel('velocity');
% 
%     subplot(3,1,3)
%     plot(accel_EE(:,3));
%     grid on
%     xlabel('time(sec)');
%     ylabel('acceleration');
% 
% end

%% Collision Check
% Create two platforms
platform1 = collisionBox(0.5,0.5,0.15);
platform1.Pose = trvec2tform([-0.5 0.4 0]);

platform2 = collisionBox(0.5,0.5,0.15);
platform2.Pose = trvec2tform([0.5 0.2 0]);

% Add a light fixture, modeled as a sphere
lightFixture1 = collisionSphere(0.05);
lightFixture1.Pose = trvec2tform([0.1 -0.65 0.7]);

lightFixture2 = collisionSphere(0.05);
lightFixture2.Pose = trvec2tform([0.2 0 0.7]);

% Store in a cell array for collision-checking
worldCollisionArray = {lightFixture1};

% Visualize the environment using a helper function that iterates through the collision array.
ax = exampleHelperVisualizeCollisionEnvironment(worldCollisionArray);

% Initialize outputs
inCollision = false(size(pos_EE,1), 1); % Check whether each pose is in collision
worldCollisionPairIdx = cell(size(pos_EE,1),1); % Provide the bodies that are in collision

%% robot visualization
panda = loadrobot("frankaEmikaPanda", "DataFormat", "row");
removeBody(panda, "panda_rightfinger");
removeBody(panda, "panda_leftfinger");
removeBody(panda, "panda_hand");

ik = inverseKinematics('RigidBodyTree', panda, 'SolverAlgorithm', 'BFGSGradientProjection');
weights = [1,1,1,1,1,1];
homeguess = panda.homeConfiguration;
initialguess = [0,0,0,-pi/2,0,pi/3,-pi/2];
randomguess = randomConfiguration(panda);

MtoH_R1 = [0, 1, 0, pi];
MtoH_R2 = [0, 0, 1, pi/4];
for i = 1:size(pos_EE,1)
    if i>1
        pos_joint(i,:) = ik('panda_link8',trvec2tform(pos_EE(i,:))*axang2tform(MtoH_R1(1,:))*axang2tform(MtoH_R2(1,:)),weights,pos_joint(i-1,:));
    else
        pos_joint(i,:) = ik('panda_link8',trvec2tform(pos_EE(i,:))*axang2tform(MtoH_R1(1,:))*axang2tform(MtoH_R2(1,:)),weights,initialguess);
    end
    lightFixture1.Pose(2,4) = lightFixture1.Pose(2,4) + 0.015;
    worldCollisionArray = {lightFixture1};
    show(worldCollisionArray{1});

    % Obstacle Avoidance
    [inCollision(i),sepDist] = checkCollision(panda,pos_joint(i,:),worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");
    [bodyIdx,worldCollisionObjIdx] = find(isnan(sepDist)); % Find collision pairs
    worldCollidingPairs = [bodyIdx,worldCollisionObjIdx]; 
    worldCollisionPairIdx{i} = worldCollidingPairs;
    
    show(panda, pos_joint(i,:),'PreservePlot',false,'visuals','on','collision','off');
    drawnow;
    pause(tf/size(pos_EE,1)); 
end

%%%%% check collision %%%%%
isTrajectoryInCollision = any(inCollision);
collidingIdx1 = find(inCollision,1);
collidingIdx2 = find(inCollision,1,"last");

if isempty(collidingIdx1)==0 && isempty(collidingIdx2)==0
    % Identify the colliding rigid bodies.
    collidingBodies1 = worldCollisionPairIdx{collidingIdx1}*[1 0]';
    collidingBodies2 = worldCollisionPairIdx{collidingIdx2}*[1 0]';
    
    % Visualize the environment.
    ax = exampleHelperVisualizeCollisionEnvironment(worldCollisionArray);
    
    % Add the robotconfigurations & highlight the colliding bodies.
    figure(2)
    show(panda,pos_joint(collidingIdx1, :),"Parent",ax,"PreservePlot",false);
    exampleHelperHighlightCollisionBodies(panda,collidingBodies1 + 1,ax);
    show(panda,pos_joint(collidingIdx2, :),"Parent"',ax);
    exampleHelperHighlightCollisionBodies(panda,collidingBodies2 + 1,ax);
else
    disp("There is not any collision")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Visualization of EE
x = linspace(0,tf,size(pos_EE,1));
y1_range=[-2.8973, 2.8973];
y2_range=[-1.7628, 1.7628];
y3_range=[-2.8973, 2.8973];
y4_range=[-3.0718, -0.0698];
y5_range=[-2.8973, 2.8973];
y6_range=[-0.0175, 3.7525];
y7_range=[-2.8973, 2.8973];

figure(3)
subplot(7,1,1)
plot(x, pos_joint(:,1));
grid on
hold on
yline(y1_range, '--r', {'Min', 'Max'})
ylabel('position-1');

subplot(7,1,2)
plot(x, pos_joint(:,2));
grid on
hold on
yline(y2_range, '--r', {'Min', 'Max'})
ylabel('position-2');

subplot(7,1,3)
plot(x, pos_joint(:,3));
grid on
hold on
yline(y3_range, '--r', {'Min', 'Max'})
ylabel('position-3');

subplot(7,1,4)
plot(x, pos_joint(:,4));
grid on
hold on
yline(y4_range, '--r', {'Min', 'Max'})
ylabel('position-4');

subplot(7,1,5)
plot(x, pos_joint(:,5));
grid on
hold on
yline(y5_range, '--r', {'Min', 'Max'})
ylabel('position-5');

subplot(7,1,6)
plot(x, pos_joint(:,6));
grid on
hold on
yline(y6_range, '--r', {'Min', 'Max'})
ylabel('position-6');

subplot(7,1,7)
plot(x, pos_joint(:,7));
grid on
hold on
yline(y7_range, '--r', {'Min', 'Max'})
xlabel('time(sec)');
ylabel('position-7');


function [ax, figHandle] = exampleHelperVisualizeCollisionEnvironment(collisionObjectArray)
    %exampleHelperVisualizeCollisionEnvironment Visualize a set of collision objects
    %   given a cell array of collision objects, create a figure that plots
    %   these and with hold on. Return the handle to the figure's axis.
    figHandle = figure;

    % Show the first object
    show(collisionObjectArray{1});
    
    % Get axis properties and set hold
    ax = gca;
    hold all;
    
    % Show remaining objects
    for i = 2:numel(collisionObjectArray)
        show(collisionObjectArray{i}, "Parent", ax);
    end
    
    % Set axis properties
    axis equal;

end

function exampleHelperHighlightCollisionBodies(robot, collisionBodyIdx, ax)
    % This function is for internal use only and may be removed in a future release.
    
    %exampleHelperURDFCollisions Check for collisions
    %   Highlight the bodies with indices given in the COLLISIONBODYIDX matrix,
    %   where the indices start correspond to the following ordering:
    %   [robot.Base robot.Bodies]. The visualization occurs in the axes
    %   specified by AX, which must already contain a visualization of the
    %   associated rigidbodytree, given by ROBOT.
    
    %   Copyright 2019 The MathWorks, Inc.
    
    % The rigid bodies actually start at the first link, not the
    % base, so the indices have to be shifted down be 1
    validateattributes(collisionBodyIdx, {'double'}, {}, 'showCollision', 'collisionBodyIdx');
    if ~isempty(collisionBodyIdx)
        rigidBodyIdx = collisionBodyIdx-1;
    else
        rigidBodyIdx = collisionBodyIdx;
    end
    
    highlightColor = [1 0.8 0]; % Yellow
    
    for i = 1:numel(rigidBodyIdx)
        if rigidBodyIdx(i) < 0
            % Body is the base
            p = findall(ax, 'type', 'patch', 'displayname', [robot.Base.Name '_mesh']);
        else
            % Any other body
            p = findall(ax, 'type', 'patch', 'displayname', [robot.Bodies{rigidBodyIdx(i)}.Name '_mesh']);
        end
        
        if isempty(p)
            continue
        else
            p(1).FaceColor = highlightColor;
        end
    end
end


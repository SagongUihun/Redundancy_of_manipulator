global traj_num
traj_num = size(pos_EE,1);

%% Init setting of manipulator
global panda
panda=loadrobot("frankaEmikaPanda", "DataFormat", "row");
removeBody(panda, "panda_rightfinger");
removeBody(panda, "panda_leftfinger");
removeBody(panda, "panda_hand");
removeBody(panda, "panda_link8");

weights = [1,1,1,1,1,1];

homeguess = panda.homeConfiguration;
initialguess = [0,0,0,-pi/2,0,pi/3,-pi/2];
randomguess = randomConfiguration(panda);

%% Collision Check
% Create two platforms
platform1 = collisionBox(0.5,0.5,0.15);
platform1.Pose = trvec2tform([-0.5 0.4 0]);

platform2 = collisionBox(0.5,0.5,0.15);
platform2.Pose = trvec2tform([0.5 0.2 0]);

% Add a light fixture, modeled as a sphere
lightFixture1 = collisionSphere(0.1);
lightFixture1.Pose = trvec2tform([-0.3 -0.3 0.4]);

% Store in a cell array for collision-checking
worldCollisionArray = {lightFixture1};

% Visualize the environment using a helper function that iterates through the collision array.
ax = exampleHelperVisualizeCollisionEnvironment(worldCollisionArray);

% Initialize outputs
inCollision = false(size(pos_EE,1), 1); % Check whether each pose is in collision
worldCollisionPairIdx = cell(size(pos_EE,1),1); % Provide the bodies that are in collision

%% WLN
joint_limit(1:7,1:2) = [-2.8973, 2.8973;
                        -1.7628, 1.7628;
                        -2.8973, 2.8973;
                        -3.0718, -0.0698;
                        -2.8973, 2.8973;
                        -0.0175, 3.7525;
                        -2.8973, 2.8973;];

rotm =[-1,0,0;
        0,1,0;
        0,0,-1];
rotXYZ=rotm2eul(rotm,'XYZ');

%%%%%%%%%%%% init_step %%%%%%%%%%%%
x(1:3,1)=pos_EE(1,:)';
x(4:6,1)=rotXYZ';
x_dot(1:3,1)=vel_EE(1,:)';
x_dot(4:6,1)=[0;0;0];

WLN_q(1:7,1)=pos_joint(1,:)';
% WLN_q(:,1)=check_over_joint_limit(WLN_q(:,1), joint_limit);

J = get_jacobian(WLN_q(:,1));
J_pinv = pinv(J,1e-4);

% Obstacle Avoidance
[inCollision(1),sepDist] = checkCollision(panda,WLN_q(:,1)',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");
[bodyIdx,worldCollisionObjIdx] = find(isnan(sepDist)); % Find collision pairs
worldCollidingPairs = [bodyIdx,worldCollisionObjIdx]; 
worldCollisionPairIdx{1} = worldCollidingPairs;

H_OA = performance_criterion_OA(sepDist);
init_g_H_OA = [0 0 0 0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
W=eye(7,7);
k=0.15; % to be maximized

% joint velocity 계산 
gen_J_pinv=inv(W)*transpose(J)*inv(J*inv(W)*transpose(J));
WLN_q_dot(:,1) = gen_J_pinv*x_dot + k*(eye(7,7) - gen_J_pinv*J)*transpose(init_g_H_OA);

now_g_H_OA = init_g_H_OA;
%%%%%%%%%%%% next_step %%%%%%%%%%%%
for i=2:size(pos_joint,1)
    x(1:3,i)=pos_EE(i,:)';
    x(4:6,i)=rotXYZ';
    x_dot(1:3,1)=vel_EE(i,:)';
    x_dot(4:6,1)=[0;0;0];

    J = get_jacobian(WLN_q(:,i-1));
    J_pinv = pinv(J,1e-4);
    
    % update generalized_pseudo_inverse of Jacobian
    gen_J_pinv=inv(W)*transpose(J)*inv(J*inv(W)*transpose(J));

    % WLN definition
    WLN_q_dot(:,i) = gen_J_pinv*x_dot + k*(eye(7,7) - gen_J_pinv*J)*transpose(now_g_H_OA);

    % WLN_q 만들어주기
    WLN_q(:,i)=WLN_q(:,i-1)+WLN_q_dot(:,i)*time_step;
%     WLN_q(:,i)=check_over_joint_limit(WLN_q(:,i), joint_limit);

    % Update H, g_H
    H_OA = performance_criterion_OA(sepDist);
    now_g_H_OA = get_gradient_H_OA(worldCollisionArray, WLN_q(:,i-1), WLN_q(:,i));

    %%%%% check collision %%%%%
    [inCollision(i),sepDist] = checkCollision(panda,WLN_q(:,i)',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");
    
    [bodyIdx,worldCollisionObjIdx] = find(isnan(sepDist)); % Find collision pairs
    worldCollidingPairs = [bodyIdx,worldCollisionObjIdx]; 
    worldCollisionPairIdx{i} = worldCollidingPairs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% visualization
% figure(6)
for i=1:size(WLN_q,2)
    show(panda, WLN_q(:,i)','PreservePlot',false,'visuals','on','collision','off');
    drawnow;
    pause(tf/size(WLN_q,2)); 
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
    show(panda,WLN_q(:,collidingIdx1)',"Parent",ax,"PreservePlot",false);
    exampleHelperHighlightCollisionBodies(panda,collidingBodies1 + 1,ax);
    show(panda,WLN_q(:,collidingIdx2)',"Parent"',ax);
    exampleHelperHighlightCollisionBodies(panda,collidingBodies2 + 1,ax);
else
    disp("There is not any collision")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WLN_q
x = linspace(0,tf,size(pos_EE,1));
y1_range=[-2.8973, 2.8973];
y2_range=[-1.7628, 1.7628];
y3_range=[-2.8973, 2.8973];
y4_range=[-3.0718, -0.0698];
y5_range=[-2.8973, 2.8973];
y6_range=[-0.0175, 3.7525];
y7_range=[-2.8973, 2.8973];

figure(7)
subplot(7,1,1)
plot(x, WLN_q(1,:));
grid on
hold on
yline(y1_range, '--r', {'Min', 'Max'})
ylabel('position-1');

subplot(7,1,2)
plot(x, WLN_q(2,:));
grid on
hold on
yline(y2_range, '--r', {'Min', 'Max'})
ylabel('position-2');

subplot(7,1,3)
plot(x, WLN_q(3,:));
grid on
hold on
yline(y3_range, '--r', {'Min', 'Max'})
ylabel('position-3');

subplot(7,1,4)
plot(x, WLN_q(4,:));
grid on
hold on
yline(y4_range, '--r', {'Min', 'Max'})
ylabel('position-4');

subplot(7,1,5)
plot(x, WLN_q(5,:));
grid on
hold on
yline(y5_range, '--r', {'Min', 'Max'})
ylabel('position-5');

subplot(7,1,6)
plot(x, WLN_q(6,:));
grid on
hold on
yline(y6_range, '--r', {'Min', 'Max'})
ylabel('position-6');

subplot(7,1,7)
plot(x, WLN_q(7,:));
grid on
hold on
yline(y7_range, '--r', {'Min', 'Max'})
xlabel('time(sec)');
ylabel('position-7');

%% WLN_q_dot
y1_range=[-2.175, 2.175];
y2_range=[-2.175, 2.175];
y3_range=[-2.175, 2.175];
y4_range=[-2.175, 2.175];
y5_range=[-2.610, 2.610];
y6_range=[-2.610, 2.610];
y7_range=[-2.610, 2.610];

figure(8)
subplot(7,1,1)
plot(x, WLN_q_dot(1,:));
grid on
hold on
yline(y1_range, '--r', {'Min', 'Max'})
ylabel('position-1');

subplot(7,1,2)
plot(x, WLN_q_dot(2,:));
grid on
hold on
yline(y2_range, '--r', {'Min', 'Max'})
ylabel('position-2');

subplot(7,1,3)
plot(x, WLN_q_dot(3,:));
grid on
hold on
yline(y3_range, '--r', {'Min', 'Max'})
ylabel('position-3');

subplot(7,1,4)
plot(x, WLN_q_dot(4,:));
grid on
hold on
yline(y4_range, '--r', {'Min', 'Max'})
ylabel('position-4');

subplot(7,1,5)
plot(x, WLN_q_dot(5,:));
grid on
hold on
yline(y5_range, '--r', {'Min', 'Max'})
ylabel('position-5');

subplot(7,1,6)
plot(x, WLN_q_dot(6,:));
grid on
hold on
yline(y6_range, '--r', {'Min', 'Max'})
ylabel('position-6');

subplot(7,1,7)
plot(x, WLN_q_dot(7,:));
grid on
hold on
yline(y7_range, '--r', {'Min', 'Max'})
xlabel('time(sec)');
ylabel('position-7');


function TF = TF_matrix(i, DH)
    d = DH(i,1);
    theta = DH(i,2);
    a = DH(i,3);
    alpha = DH(i,4);
    
    TF = [ cos(theta), -cos(alpha)*sin(theta),  sin(alpha)*sin(theta), a*cos(theta) ;
           sin(theta),  cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta) ;
                0,            sin(alpha),              cos(alpha),          d       ;
                0,                0,                       0,               1      ];       
end
function Jacob = jacobian(T76, T65, T54, T43, T32, T21, T10)
    T75 = T76*T65;
    T74 = T75*T54;
    T73 = T74*T43;
    T72 = T73*T32;
    T71 = T72*T21;
    T70 = T71*T10;
    
    o0 = [0; 0; 0];
    o1 = T76(1:3,4);
    o2 = T75(1:3,4);
    o3 = T74(1:3,4);
    o4 = T73(1:3,4);
    o5 = T72(1:3,4);
    o6 = T71(1:3,4);
    o7 = T70(1:3,4);
    
    z0 = [0; 0; 1];
    z1 = T76(1:3,3);
    z2 = T75(1:3,3);
    z3 = T74(1:3,3);
    z4 = T73(1:3,3);
    z5 = T72(1:3,3);
    z6 = T71(1:3,3);
    z7 = T70(1:3,3);
    
    Jacob(1:3, 1) = cross(z0, o7 - o0);
    Jacob(1:3, 2) = cross(z1, o7 - o1);
    Jacob(1:3, 3) = cross(z2, o7 - o2);
    Jacob(1:3, 4) = cross(z3, o7 - o3);
    Jacob(1:3, 5) = cross(z4, o7 - o4);
    Jacob(1:3, 6) = cross(z5, o7 - o5);
    Jacob(1:3, 7) = cross(z6, o7 - o6);
    
    Jacob(4:6, 1) = z0;
    Jacob(4:6, 2) = z1;
    Jacob(4:6, 3) = z2;
    Jacob(4:6, 4) = z3;
    Jacob(4:6, 5) = z4;
    Jacob(4:6, 6) = z5;
    Jacob(4:6, 7) = z6;
end
function Jacob07 = get_jacobian(config_q)
    theta1 = config_q(1,1);
    theta2 = config_q(2,1);
    theta3 = config_q(3,1);
    theta4 = config_q(4,1);
    theta5 = config_q(5,1);
    theta6 = config_q(6,1);
    theta7 = config_q(7,1);

    a = 0.333;
    b = 0.316;
    c = 0.384;
    d = 0.088;
    e = 0.107;
    f = 0.0825;

    DH =  [ a,  theta1,    0,  -pi/2;
            0,  theta2,    0,   pi/2;
            b,  theta3,    f,   pi/2;
            0,  theta4,    -f, -pi/2;
            c,  theta5,    0,   pi/2;
            0,  theta6,    d,   pi/2;
            e,  theta7 + pi/4,    0,   0  ];

    T01 = TF_matrix(1, DH);
    T12 = TF_matrix(2, DH);
    T23 = TF_matrix(3, DH);
    T34 = TF_matrix(4, DH);
    T45 = TF_matrix(5, DH);
    T56 = TF_matrix(6, DH);
    T67 = TF_matrix(7, DH);

    Jacob07 = jacobian(T01, T12, T23, T34, T45, T56, T67);
end
function H = performance_criterion_JL(joint_limit, config_q)
    H1=(joint_limit(1,2)-joint_limit(1,1))*(joint_limit(1,2)-joint_limit(1,1))/(4*(joint_limit(1,2)-config_q(1,1))*(config_q(1,1)-joint_limit(1,1)));
    H2=(joint_limit(2,2)-joint_limit(2,1))*(joint_limit(2,2)-joint_limit(2,1))/(4*(joint_limit(2,2)-config_q(2,1))*(config_q(2,1)-joint_limit(2,1)));
    H3=(joint_limit(3,2)-joint_limit(3,1))*(joint_limit(3,2)-joint_limit(3,1))/(4*(joint_limit(3,2)-config_q(3,1))*(config_q(3,1)-joint_limit(3,1)));
    H4=(joint_limit(4,2)-joint_limit(4,1))*(joint_limit(4,2)-joint_limit(4,1))/(4*(joint_limit(4,2)-config_q(4,1))*(config_q(4,1)-joint_limit(4,1)));
    H5=(joint_limit(5,2)-joint_limit(5,1))*(joint_limit(5,2)-joint_limit(5,1))/(4*(joint_limit(5,2)-config_q(5,1))*(config_q(5,1)-joint_limit(5,1)));
    H6=(joint_limit(6,2)-joint_limit(6,1))*(joint_limit(6,2)-joint_limit(6,1))/(4*(joint_limit(6,2)-config_q(6,1))*(config_q(6,1)-joint_limit(6,1)));
    H7=(joint_limit(7,2)-joint_limit(7,1))*(joint_limit(7,2)-joint_limit(7,1))/(4*(joint_limit(7,2)-config_q(7,1))*(config_q(7,1)-joint_limit(7,1)));

    H = H1 + H2 + H3 + H4 + H5 + H6 + H7;
end
function g_H= get_gradient_H_JL(joint_limit, config_q)
    g_H1 = ((joint_limit(1,2)-joint_limit(1,1))*(joint_limit(1,2)-joint_limit(1,1))*(2*config_q(1,1)-joint_limit(1,2)-joint_limit(1,1)))/(4*(joint_limit(1,2)-config_q(1,1))*(joint_limit(1,2)-config_q(1,1))*(config_q(1,1)-joint_limit(1,1))*(config_q(1,1)-joint_limit(1,1)));
    g_H2 = ((joint_limit(2,2)-joint_limit(2,1))*(joint_limit(2,2)-joint_limit(2,1))*(2*config_q(2,1)-joint_limit(2,2)-joint_limit(2,1)))/(4*(joint_limit(2,2)-config_q(2,1))*(joint_limit(2,2)-config_q(2,1))*(config_q(2,1)-joint_limit(2,1))*(config_q(2,1)-joint_limit(2,1)));
    g_H3 = ((joint_limit(3,2)-joint_limit(3,1))*(joint_limit(3,2)-joint_limit(3,1))*(2*config_q(3,1)-joint_limit(3,2)-joint_limit(3,1)))/(4*(joint_limit(3,2)-config_q(3,1))*(joint_limit(3,2)-config_q(3,1))*(config_q(3,1)-joint_limit(3,1))*(config_q(3,1)-joint_limit(3,1)));
    g_H4 = ((joint_limit(4,2)-joint_limit(4,1))*(joint_limit(4,2)-joint_limit(4,1))*(2*config_q(4,1)-joint_limit(4,2)-joint_limit(4,1)))/(4*(joint_limit(4,2)-config_q(4,1))*(joint_limit(4,2)-config_q(4,1))*(config_q(4,1)-joint_limit(4,1))*(config_q(4,1)-joint_limit(4,1)));
    g_H5 = ((joint_limit(5,2)-joint_limit(5,1))*(joint_limit(5,2)-joint_limit(5,1))*(2*config_q(5,1)-joint_limit(5,2)-joint_limit(5,1)))/(4*(joint_limit(5,2)-config_q(5,1))*(joint_limit(5,2)-config_q(5,1))*(config_q(5,1)-joint_limit(5,1))*(config_q(5,1)-joint_limit(5,1)));
    g_H6 = ((joint_limit(6,2)-joint_limit(6,1))*(joint_limit(6,2)-joint_limit(6,1))*(2*config_q(6,1)-joint_limit(6,2)-joint_limit(6,1)))/(4*(joint_limit(6,2)-config_q(6,1))*(joint_limit(6,2)-config_q(6,1))*(config_q(6,1)-joint_limit(6,1))*(config_q(6,1)-joint_limit(6,1)));
    g_H7 = ((joint_limit(7,2)-joint_limit(7,1))*(joint_limit(7,2)-joint_limit(7,1))*(2*config_q(7,1)-joint_limit(7,2)-joint_limit(7,1)))/(4*(joint_limit(7,2)-config_q(7,1))*(joint_limit(7,2)-config_q(7,1))*(config_q(7,1)-joint_limit(7,1))*(config_q(7,1)-joint_limit(7,1)));
    
    g_H = [g_H1, g_H2, g_H3, g_H4, g_H5, g_H6, g_H7];
end

function H = performance_criterion_OA(sepDist)
    Dist=sepDist(1:7,:);
    min_d = min(Dist);
    
    k = 0.5;
    for i=1:size(Dist,1)
        b(i,1) = k^(Dist(i,1) - min_d);
        h(i,1) = b(i,1)*Dist(i,1);
    end
    
    H = sum(h);
end

function g_H = get_gradient_H_OA(worldCollisionArray, past_q, now_q)
    global traj_num
    global panda
    inCollision_past = false(traj_num, 1); 
    [inCollision_past(1),sepDist_past] = checkCollision(panda,past_q',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");
    
    % "new_q" that is added each q_i to past_q
    past_q_delta_q1 = past_q;
    past_q_delta_q1(1,1) = now_q(1,1);
    inCollision_H1 = false(traj_num, 1); 
    [inCollision_H1(1),sepDist_H1] = checkCollision(panda,past_q_delta_q1',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");
    
    past_q_delta_q2 = past_q;
    past_q_delta_q2(2,1) = now_q(2,1);
    inCollision_H2 = false(traj_num, 1); 
    [inCollision_H2(1),sepDist_H2] = checkCollision(panda,past_q_delta_q2',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");

    past_q_delta_q3 = past_q;
    past_q_delta_q3(3,1) = now_q(3,1);
    inCollision_H3 = false(traj_num, 1); 
    [inCollision_H3(1),sepDist_H3] = checkCollision(panda,past_q_delta_q3',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");

    past_q_delta_q4 = past_q;
    past_q_delta_q4(4,1) = now_q(4,1);
    inCollision_H4 = false(traj_num, 1); 
    [inCollision_H4(1),sepDist_H4] = checkCollision(panda,past_q_delta_q4',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");

    past_q_delta_q5 = past_q;
    past_q_delta_q5(5,1) = now_q(5,1);
    inCollision_H5 = false(traj_num, 1); 
    [inCollision_H5(1),sepDist_H5] = checkCollision(panda,past_q_delta_q5',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");

    past_q_delta_q6 = past_q;
    past_q_delta_q6(6,1) = now_q(6,1);
    inCollision_H6 = false(traj_num, 1); 
    [inCollision_H6(1),sepDist_H6] = checkCollision(panda,past_q_delta_q6',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");

    past_q_delta_q7 = past_q;
    past_q_delta_q7(7,1) = now_q(7,1);
    inCollision_H7 = false(traj_num, 1); 
    [inCollision_H7(1),sepDist_H7] = checkCollision(panda,past_q_delta_q7',worldCollisionArray,"IgnoreSelfCollision","on","Exhaustive","on","SkippedSelfCollisions","parent");

    % Calculate gradient_H
    g_H1 = ( performance_criterion_OA(sepDist_H1) - performance_criterion_OA(sepDist_past) ) / (now_q(1,1) - past_q(1,1));
    g_H2 = ( performance_criterion_OA(sepDist_H2) - performance_criterion_OA(sepDist_past) ) / (now_q(2,1) - past_q(2,1));
    g_H3 = ( performance_criterion_OA(sepDist_H3) - performance_criterion_OA(sepDist_past) ) / (now_q(3,1) - past_q(3,1));
    g_H4 = ( performance_criterion_OA(sepDist_H4) - performance_criterion_OA(sepDist_past) ) / (now_q(4,1) - past_q(4,1));
    g_H5 = ( performance_criterion_OA(sepDist_H5) - performance_criterion_OA(sepDist_past) ) / (now_q(5,1) - past_q(5,1));
    g_H6 = ( performance_criterion_OA(sepDist_H6) - performance_criterion_OA(sepDist_past) ) / (now_q(6,1) - past_q(6,1));
    g_H7 = ( performance_criterion_OA(sepDist_H7) - performance_criterion_OA(sepDist_past) ) / (now_q(7,1) - past_q(7,1));
    
    g_H = [g_H1, g_H2, g_H3, g_H4, g_H5, g_H6, g_H7];
end

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

function refined_q = check_over_joint_limit(config_q, joint_limit)
    refined_q = config_q;
    joint_limit(1:7,1:2) = [-2.8973, 2.8973;
                            -1.7628, 1.7628;
                            -2.8973, 2.8973;
                            -3.0718, -0.0698;
                            -2.8973, 2.8973;
                            -0.0175, 3.7525;
                            -2.8973, 2.8973;];
    if config_q(1,1) < joint_limit(1,1)
        refined_q(1,1) = joint_limit(1,1);
    elseif config_q(1,1) > joint_limit(1,2)
        refined_q(1,1) = joint_limit(1,2);
    end
    
    if config_q(2,1) < joint_limit(2,1)
        refined_q(2,1) = joint_limit(2,1);
    elseif config_q(2,1) > joint_limit(2,2)
        refined_q(2,1) = joint_limit(2,2);
    end

    if config_q(3,1) < joint_limit(3,1)
        refined_q(3,1) = joint_limit(3,1);
    elseif config_q(3,1) > joint_limit(3,2)
        refined_q(3,1) = joint_limit(3,2);
    end

    if config_q(4,1) < joint_limit(4,1)
        refined_q(4,1) = joint_limit(4,1);
    elseif config_q(4,1) > joint_limit(4,2)
        refined_q(4,1) = joint_limit(4,2);
    end

    if config_q(5,1) < joint_limit(5,1)
        refined_q(5,1) = joint_limit(5,1);
    elseif config_q(5,1) > joint_limit(5,2)
        refined_q(5,1) = joint_limit(5,2);
    end

    if config_q(6,1) < joint_limit(6,1)
        refined_q(6,1) = joint_limit(6,1);
    elseif config_q(6,1) > joint_limit(6,2)
        refined_q(6,1) = joint_limit(6,2);
    end

    if config_q(7,1) < joint_limit(7,1)
        refined_q(7,1) = joint_limit(7,1);
    elseif config_q(7,1) > joint_limit(7,2)
        refined_q(7,1) = joint_limit(7,2);
    end
end

import pygame
from pygame import mixer
import sys
import random
import os
import time

# Region "Global variables and settings"
os.chdir(sys.path[0])

# re: grid
screen_height = screen_width = 480
gridsize = screen_height/24
grid_width = screen_width/gridsize
grid_height = screen_height/gridsize

# re: moves
up = (0,-1)
down = (0,1)
left = (-1,0)
right = (1,0)

# choose colours (R,G,B)
snake_color1 = (15, 56, 15)
snake_color2 = (15, 56, 15)
food_color = (15, 56, 15)
invert_light = (20, 40, 20)

# board_light = (139, 172, 15)
board_light = (120, 155, 11)
board_dark = (130, 160, 15)

# game options
snake_length = 5
food_target = 8
result_display_time = 0.8
flash = 6

class Snake():
    def __init__(self):
        # player 1
        self.length_p1 = snake_length
        self.positions_p1 = [((60), (screen_height/2))]
        self.direction_p1 = down
        self.color_p1 = snake_color1
        self.score_p1 = 0
        self.overall_p1 = 0
        self.collision_p1 = False
        # player 2
        self.length_p2 = snake_length
        self.positions_p2 = [((screen_width-(60)), (screen_height/2))]
        self.direction_p2 = up
        self.color_p2 = snake_color2
        self.score_p2 = 0
        self.overall_p2 = 0
        self.collision_p2 = False

    def get_head_position_p1(self):
        # return the position of the head of the snake
        return self.positions_p1[0]

    def get_head_position_p2(self):
        # return the position of player 2 snake
        return self.positions_p2[0]

    def turn_p1(self, point):
        # if the snake is one block, it can move in any direction
        if self.length_p1 > 1 and (point[0]*-1, point[1]*-1) == self.direction_p1:
            return
        else:
            self.direction_p1 = point

    def turn_p2(self, point):
        # if the snake is one block, it can move in any direction
        if self.length_p2 > 1 and (point[0]*-1, point[1]*-1) == self.direction_p2:
            return
        else:
            self.direction_p2 = point

    def move_p1(self):
        # function to move the snake
        cur = self.get_head_position_p1()
        x,y = self.direction_p1
        new = (((cur[0]+(x*gridsize))%screen_width), (cur[1]+(y*gridsize))%screen_height)
        # if the game is lost
        if len(self.positions_p1) > 2 and new in self.positions_p1[2:]:
            self.collision_p1 = True
        else:
            # otherwise, move the snake
            self.positions_p1.insert(0,new)
            if len(self.positions_p1) > self.length_p1:
                self.positions_p1.pop()

    def move_p2(self):
        # function to move the snake
        cur = self.get_head_position_p2()
        x,y = self.direction_p2
        new = (((cur[0]+(x*gridsize))%screen_width), (cur[1]+(y*gridsize))%screen_height)
        # if the game is lost
        if len(self.positions_p2) > 2 and new in self.positions_p2[2:]:
            self.collision_p2 = True
        else:
            # otherwise, grow the snake
            self.positions_p2.insert(0,new)
            if len(self.positions_p2) > self.length_p2:
                self.positions_p2.pop()

    def reset(self):
        # when the game ends
        # player 1
        self.length_p1 = snake_length
        self.positions_p1 = [((60), (screen_height/2))]
        self.direction_p1 = down
        self.score_p1 = 0
        # player 2
        self.length_p2 = snake_length
        self.positions_p2 = [((screen_width-(60)), (screen_height/2))]
        self.direction_p2 = up
        self.score_p2 = 0

    def soft_reset(self):
        # when players clash heads
        # player 1
        self.positions_p1 = [((60), (screen_height/2))]
        self.direction_p1 = down
        # player 2
        self.positions_p2 = [((screen_width-(60)), (screen_height/2))]
        self.direction_p2 = up

    def draw_p1(self,surface):
        # create a draw function
        for p in self.positions_p1:
            r = pygame.Rect((p[0], p[1]), (gridsize,gridsize))
            pygame.draw.rect(surface, self.color_p1, r)
            pygame.draw.rect(surface, board_light, r, 1)

    def draw_p2(self,surface):
        # create a draw function
        for p in self.positions_p2:
            r = pygame.Rect((p[0], p[1]), (gridsize,gridsize))
            pygame.draw.rect(surface, self.color_p2, r)
            pygame.draw.rect(surface, board_light, r, 1)

    def handle_keys(self):
        # handle user input
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_UP:
                    self.turn_p2(up)
                elif event.key == pygame.K_DOWN:
                    self.turn_p2(down)
                elif event.key == pygame.K_LEFT:
                    self.turn_p2(left)
                elif event.key == pygame.K_RIGHT:
                    self.turn_p2(right)
                elif event.key == pygame.K_w:
                    self.turn_p1(up)
                elif event.key == pygame.K_a:
                    self.turn_p1(left)
                elif event.key == pygame.K_d:
                    self.turn_p1(right)
                elif event.key == pygame.K_s:
                    self.turn_p1(down)

class Food():
    def __init__(self):
        # initialise food properties
        self.position = (0,0)
        self.color = food_color

    def randomize_position(self, snake):
        # randomise position
        self.position = (random.randint(0, grid_width-1)*gridsize, random.randint(0, grid_height-1)*gridsize)
        while self.position in snake.positions_p1 or self.position in snake.positions_p2:
            self.position = (random.randint(0, grid_width-1)*gridsize, random.randint(0, grid_height-1)*gridsize)

    def draw(self, surface):
        # draw the food onto the surface
        r = pygame.Rect((self.position[0], self.position[1]), (gridsize, gridsize))
        pygame.draw.rect(surface, self.color, r)
        pygame.draw.rect(surface, board_light, r, 1)

class Game():
    def __init__(self):
        # set pygame initialisation
        pygame.mixer.init()
        pygame.mixer.music.set_volume(0.2)
        pygame.init()
        # annotate screen
        pygame.display.set_caption('Snake')
        a = pygame.image.load('assets/icon/icon.jpg')
        pygame.display.set_icon(a)
        # soundtracks
        mixer.music.load('assets/music/soundtrack/st1.wav')
        mixer.music.play(-1)
        
        # set pygame variables
        self.clock = pygame.time.Clock()
        self.screen = pygame.display.set_mode((screen_width, screen_height), 0, 32)
        # create surface
        self.surface = pygame.Surface(self.screen.get_size())
        self.surface = self.surface.convert()
        drawGrid(self.surface, False)
        # pygame fonts
        self.myfont = pygame.font.Font("assets/fonts/8-BitMadness.ttf",26, bold=True)
        self.mylargefont = pygame.font.Font("assets/fonts/8-BitMadness.ttf",36, bold=True)

def drawGrid(surface, inverted):
    if inverted:
        # function that draws the grid
        for y in range(0, int(grid_height)):
            for x in range(0, int(grid_width)):
                if (x+y)%2 == 0:
                    r = pygame.Rect((x*gridsize, y*gridsize), (gridsize,gridsize))
                    pygame.draw.rect(surface, invert_light, r)
                else:
                    rr = pygame.Rect((x*gridsize, y*gridsize), (gridsize,gridsize))
                    pygame.draw.rect(surface, snake_color1, rr)
    else:
        # function that draws the grid
        for y in range(0, int(grid_height)):
            for x in range(0, int(grid_width)):
                if (x+y)%2 == 0:
                    r = pygame.Rect((x*gridsize, y*gridsize), (gridsize,gridsize))
                    pygame.draw.rect(surface, board_light, r)
                else:
                    rr = pygame.Rect((x*gridsize, y*gridsize), (gridsize,gridsize))
                    pygame.draw.rect(surface, board_dark, rr)

def multiplayer():
    game = Game()
    snake = Snake()
    food = Food()

    while (True):
        game.clock.tick(5 + round((snake.score_p1 + snake.score_p2)/8, 1))
        snake.handle_keys()
        drawGrid(game.surface, False)
        snake.move_p1()
        snake.move_p2()

        # snakes clash heads, no points awarded, keep playing
        if snake.get_head_position_p1() == snake.get_head_position_p2():
            snake.soft_reset()
            centre_display("No-one's point", True)
            centre_display("Keep playing!", True)
        
        # snake 1 self collision
        if snake.collision_p1:
            snake.reset()
            snake.overall_p2 += 1
            game.screen.blit(game.surface, (0,0))
            show_scores(False)
            centre_display("Player 2's point", True)
            centre_display("Player 1 got tangled up", True)
            snake.collision_p1 = False
        
        # snake 2 self collision
        if snake.collision_p2:
            snake.reset()
            snake.overall_p1 += 1
            game.screen.blit(game.surface, (0,0))
            show_scores(False)
            centre_display("Player 1's point", True)
            centre_display("Player 2 got tangled up", True)
            snake.collision_p2 = False

        # snake 1 eats a piece of food
        if snake.get_head_position_p1() == food.position:
            eating_sound = mixer.Sound("assets/music/sound_effects/eat/eat.mp3")
            eating_sound.play()
            snake.length_p1 += 1
            snake.score_p1 += 1
            food.randomize_position(snake)

        # snake 2 eats a piece of food
        if snake.get_head_position_p2() == food.position:
            eating_sound = mixer.Sound("assets/music/sound_effects/eat/eat.mp3")
            eating_sound.play()
            snake.length_p2 += 1
            snake.score_p2 += 1
            food.randomize_position(snake)

        # player 1 hits player 2
        if snake.get_head_position_p1() in snake.positions_p2[1:]:
            snake.reset()
            snake.overall_p2 += 1
            centre_display("Player 2's point", True)
            food.randomize_position(snake)
        
        # player 2 hits player 1
        if snake.get_head_position_p2() in snake.positions_p1[1:]:
            snake.reset()
            snake.overall_p1 += 1
            centre_display("Player 1's point", True)
            food.randomize_position(snake)

        # player 1 reaches the food target
        if snake.score_p1 == food_target:
            snake.reset()
            snake.overall_p1 += 1
            centre_display("Player 1's point", False)
            centre_display("Player 1 ate all the food", True)

        # player 2 reaches the food target
        if snake.score_p2 == food_target:
            snake.reset()
            snake.overall_p2 += 1
            centre_display("Player 2's point", False)
            centre_display("Player 2 ate all the food", True)

        snake.draw_p1(game.surface)
        snake.draw_p2(game.surface)
        food.draw(game.surface)

        game.screen.blit(game.surface, (0,0))

        def show_scores(inverted):
            if inverted:
                # add vs scores to screen
                score = game.mylargefont.render("(" + str(snake.score_p1) + ")" + 
                    str(snake.overall_p1) + ":" + str(snake.overall_p2) + 
                    "(" + str(snake.score_p2) + ")", 1, board_light)
                score_rect = score.get_rect(center=(screen_width/2, screen_height/25))
                game.screen.blit(score, score_rect)
            else:
                # add vs scores to screen
                score = game.mylargefont.render("(" + str(snake.score_p1) + ")" + 
                    str(snake.overall_p1) + ":" + str(snake.overall_p2) + 
                    "(" + str(snake.score_p2) + ")", 1, snake_color1)
                score_rect = score.get_rect(center=(screen_width/2, screen_height/25))
                game.screen.blit(score, score_rect)

        def centre_display(text, sound):
            if sound:
                eating_sound = mixer.Sound("assets/music/sound_effects/eat/eat.mp3")
                eating_sound.play()
            for i in range(0, 8):
                if (i % 2) == 0:
                    drawGrid(game.surface, True)
                    # blit screen
                    game.screen.blit(game.surface, (0,0))
                    # display scores
                    show_scores(True)
                    # text render
                    status = game.mylargefont.render(text, 1, board_light)
                    status_rect = status.get_rect(center=(screen_width/2, screen_height/2))
                    game.screen.blit(status, status_rect)
                    pygame.display.update()
                    # screen pause
                    time.sleep(result_display_time/flash)
                    pygame.display.update()
                else:
                    drawGrid(game.surface, False)
                    # blit screen
                    game.screen.blit(game.surface, (0,0))
                    # display scores
                    show_scores(False)
                    # text render
                    status = game.mylargefont.render(text, 1, snake_color1)
                    status_rect = status.get_rect(center=(screen_width/2, screen_height/2))
                    game.screen.blit(status, status_rect)
                    pygame.display.update()
                    # screen pause
                    time.sleep(result_display_time/flash)
                    pygame.display.update()

        # show scores on screen
        show_scores(False)
        # update the screen display
        pygame.display.update()

def singleplayer():
    game = Game()
    snake = Snake()
    food = Food()
    
    # randomise the food position
    food.randomize_position(snake)

    # load highscore data
    hstext = open("data/hstext.txt")
    highscore = int(hstext.read())
    init_highscore = highscore
        
    while(True):
        game.clock.tick(5 + round((snake.score_p1 + snake.score_p2)/8, 1))
        snake.handle_keys()
        drawGrid(game.surface, False)
        snake.move_p2()
        
        # snake self collision
        if snake.collision_p2:
            if snake.score_p2 > init_highscore:
                centre_display("New High Score!", True)
            centre_display("You scored " + str(snake.score_p2), True)
            snake.reset()
            game.screen.blit(game.surface, (0,0))
            snake.collision_p2 = False

            
        # snake eats a piece of food
        if snake.get_head_position_p2() == food.position:
            eating_sound = mixer.Sound("assets/music/sound_effects/eat/eat.mp3")
            eating_sound.play()
            snake.length_p2 += 1
            snake.score_p2 += 1
            food.randomize_position(snake)
            
            # if this is a new highscore
            if snake.score_p2 > highscore:
                highscore = snake.score_p2
                hswrite = open("data/hstext.txt", "w")
                hswrite.write(str(highscore))
                    
        # show the food and snake
        snake.draw_p2(game.surface)
        food.draw(game.surface)
        
        game.screen.blit(game.surface, (0,0))
        score = game.myfont.render("Score {0}".format(snake.score_p2), 1, snake_color2)
        game.screen.blit(score, (5,25))
        high_score = game.myfont.render("High score {0}".format(highscore), 1, snake_color2)
        game.screen.blit(high_score, (5,5))
        
        # update the screen display
        pygame.display.update()
        
        def centre_display(text, sound):
            if sound:
                eating_sound = mixer.Sound("assets/music/sound_effects/eat/eat.mp3")
                eating_sound.play()
            for i in range(0, 8):
                if (i % 2) == 0:
                    drawGrid(game.surface, True)
                    # blit screen
                    game.screen.blit(game.surface, (0,0))
                    # text render
                    status = game.mylargefont.render(text, 1, board_light)
                    status_rect = status.get_rect(center=(screen_width/2, screen_height/2))
                    game.screen.blit(status, status_rect)
                    pygame.display.update()
                    # screen pause
                    time.sleep(result_display_time/flash)
                    pygame.display.update()
                else:
                    drawGrid(game.surface, False)
                    # blit screen
                    game.screen.blit(game.surface, (0,0))
                    # text render
                    status = game.mylargefont.render(text, 1, snake_color1)
                    status_rect = status.get_rect(center=(screen_width/2, screen_height/2))
                    game.screen.blit(status, status_rect)
                    pygame.display.update()
                    # screen pause
                    time.sleep(result_display_time/flash)
                    pygame.display.update()
    
# This notation means main() will only run when this script is excecuted
if __name__ == "__main__":
    singleplayer()
    # multiplayer()

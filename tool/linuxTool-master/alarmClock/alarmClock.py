#! /usr/bin/python
# -*-coding:utf-8-*-

# !/bin/bash
# 52 05 * * * /Users/apple/Desktop/c/project/linuxTool/alarmClock/alarmClock.py > /dev/console
# 52 05 * * * /Users/apple/Desktop/c/project/linuxTool/alarmClock/alarmClock.py > /dev/null 2>&1  不输出日志
#echo "定时任务";

import os
import time


file = r'/Users/apple/Music/iTunes/五色石南叶-此处长安.mp3'     # 音乐文件地址
music_path = r'/Users/apple/Music/QQ音乐'                     # 音乐文件目录


def get_music(path):
    '''
    read all music files from the folder
    :param path:
    :return:
    '''
    row_filenames = os.listdir(path)

    music_files = []
    for file_name in row_filenames:
        if(file_name.lower().endswith('.mp3') or file_name.lower().endswith('.ogg')):
            music_files.append(os.path.join(music_path, file_name))

    return sorted(music_files)


def worker(file):
    '''
    play the music
    :param file: file path
    :return:
    '''
    import pygame
    pygame.mixer.init()  # 初始化音频部分
    if(len(file) > 0):
        for m in file:
            if not os.path.exists(m):
                print('File doesn\'t exist')
                print(pygame.ver)
                return
            # 载入将要播放的音乐文件
            track = pygame.mixer_music.load(file[2])   # 遇到一个.mp3的文件到这一步会卡死，进程在走，但是没有加载进来文件，也没有错误返回，如果有大神解决了这个问题，希望不吝告知
            while pygame.mixer.music.get_busy() == 0:  # 判断是否在播放音乐
                pygame.mixer_music.play(loops=4)  # 播放载入的音乐 loops代表播放的重复次数
                time.sleep(600)
                pygame.mixer_music.stop()  # 停止播放
                break

    else:
        print('loading music files failed')



if __name__ == '__main__':
    musics = get_music(music_path)
    worker(musics)
    # if(len(musics) > 0):
    #     for m in musics:
    #         print(m)
    #         worker(m)
    # else:
    #     worker(file)


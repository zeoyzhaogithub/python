# mac 写定时任务
## 起床闹钟自动播放音乐

如何启用crontab
首先，既然OS X的定时任务都是由launchctl进行管理，先确定是否有cron任务：

    $LaunchAgents sudo launchctl list | grep cron
    -       0 com.vix.cron

果然，在里面，接下来检查配置项

    $LaunchAgents locate com.vix.cron

    WARNING: The locate database (/var/db/locate.database) does not exist.
    To create the database, run the following command:

    sudo launchctl load -w /System/Library/LaunchDaemons/com.apple.locate.plist

    Please be aware that the database can take some time to generate; once
    the database has been created, this message will no longer appear.

第一次使用locate命令报错了，按照提示，运行命令：

    $LaunchAgents sudo launchctl load -w /System/Library/LaunchDaemons/com.apple.locate.plist

运行这个命令需要一点时间，且locate命令每次运行前都需要更新数据库

    $LaunchAgents sudo /usr/libexec/locate.updatedb

运行该命令不在报错，则数据库创建成功

重新运行命令：

    $LaunchAgents locate com.vix.cron
    /System/Library/LaunchDaemons/com.vix.cron.plist

进入文件

    $LaunchAgents cat /System/Library/LaunchDaemons/com.vix.cron.plist

    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN"
        "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
    <plist version="1.0">
    <dict>
        <key>Label</key>
        <string>com.vix.cron</string>
        <key>ProgramArguments</key>
        <array>
            <string>/usr/sbin/cron</string>
        </array>
        <key>KeepAlive</key>
        <dict>
            <key>PathState</key>
            <dict>
                <key>/etc/crontab</key>
                <true/>
            </dict>
        </dict>
        <key>QueueDirectories</key>
        <array>
            <string>/usr/lib/cron/tabs</string>
        </array>
        <key>EnableTransactions</key>
        <true/>
    </dict>
    </plist>

注意：文件里面有一个KeepAlive的条件是：/etc/crontab是否存在，如果该文件不存在，会导致定时任务不起作用

    $LaunchAgents ll /etc/crontab

1）如果提示ll command not find

    $LaunchAgents vim ~/.bash_profile

在该文件中加入下面内容

    14 # mac命令缩写别名
    15 alias ll='ls -alF'
    16 alias la='ls -A'
    17 alias l='ls -CF'

重启文件，使命令生效，

    $LaunchAgents source ~/.bash_profile

2）如果提示：

    ls: /etc/crontab: No such file or directory

文件不存在，先创建文件

    $LaunchAgents sudo touch /etc/crontab

运行脚本：

提示 Permission denied
文件没有运行权限，修改文件权限

    chmod 777 cronStudy.sh

编写脚本文件 alarmClock.py

添加定时任务

    crontab -e

进入文件，添加定时任务

    52 5 * * * /Users/apple/Desktop/c/stg/linux/cronStudy.py > /dev/null 2>&1

定时每天早晨5点52开始启动程序，关于crontab命令，请自行百度或谷歌


